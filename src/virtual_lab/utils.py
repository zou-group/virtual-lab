"""Contains useful utility functions."""

import json
import time
import urllib.parse
from pathlib import Path

import requests
import tiktoken
from openai.types.chat import ChatCompletionMessageParam
from openai.types.chat.chat_completion_message_tool_call import ChatCompletionMessageToolCall

from virtual_lab.constants import (
    DEFAULT_FINETUNING_EPOCHS,
    MODEL_TO_INPUT_PRICE_PER_TOKEN,
    MODEL_TO_OUTPUT_PRICE_PER_TOKEN,
    FINETUNING_MODEL_TO_TRAINING_PRICE_PER_TOKEN,
    PUBMED_TOOL_NAME,
    RETRACTION_TOOL_NAME,
)
from virtual_lab.prompts import format_references


def get_pubmed_central_article(pmcid: str, abstract_only: bool = False) -> tuple[str | None, list[str] | None]:
    """Gets the title and content (abstract or full text) of a PubMed Central article given a PMC ID.

    Note: This only returns main text, ignoring tables, figures, and references.

    :param pmcid: The PMC ID of the article.
    :param abstract_only: Whether to return only the abstract instead of the full text.
    :return: The title and content (abstract or full text of the article as a list of paragraphs)
        or None if the article is not found.
    """
    # Get article from PMC ID in JSON form
    text_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_JSON/PMC{pmcid}/unicode"
    response = requests.get(text_url)
    response.raise_for_status()

    # Try to parse JSON
    try:
        article = response.json()
    except json.JSONDecodeError:
        return None, None

    # Get document
    document = article[0]["documents"][0]

    # Get title
    title = next(passage["text"] for passage in document["passages"] if passage["infons"]["section_type"] == "TITLE")

    # Get relevant passages
    passages = [passage for passage in document["passages"] if passage["infons"]["type"] in {"abstract", "paragraph"}]

    # Get abstract or full text of article (excluding references)
    if abstract_only:
        passages = [passage for passage in passages if passage["infons"]["section_type"] in ["ABSTRACT"]]
    else:
        passages = [
            passage
            for passage in passages
            if passage["infons"]["section_type"] in ["ABSTRACT", "INTRO", "RESULTS", "DISCUSS", "CONCL", "METHODS"]
        ]

    # Get content
    content = [passage["text"] for passage in passages]

    return title, content


def run_pubmed_search(query: str, num_articles: int = 3, abstract_only: bool = False) -> str:
    """Runs a PubMed search, returning the full text of the top matching article.

    :param query: The query to search PubMed with.
    :param num_articles: The number of articles to search for.
    :param abstract_only: Whether to return only the abstract instead of the full text.
    :return: The full text of the top matching article.
    """
    # Print search query
    print(
        f'Searching PubMed Central for {num_articles} articles ({"abstracts" if abstract_only else "full text"}) with query: "{query}"'
    )

    # Perform PubMed Central search for query to get PMC ID
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term={urllib.parse.quote_plus(query)}&retmax={2 * num_articles}&retmode=json&sort=relevance"
    response = requests.get(search_url)
    response.raise_for_status()
    pmcids_found = response.json()["esearchresult"]["idlist"]

    # Loop through top articles
    texts = []
    titles = []
    pmcids = []

    for pmcid in pmcids_found:
        # Break if reached desired number of articles
        if len(pmcids) >= num_articles:
            break

        title, content = get_pubmed_central_article(
            pmcid=pmcid,
            abstract_only=abstract_only,
        )

        if title is None:
            continue

        texts.append(f"PMCID = {pmcid}\n\nTitle = {title}\n\n{'\n\n'.join(content or [])}")
        titles.append(title)
        pmcids.append(pmcid)

    # Print articles found
    article_count = len(texts)

    print(f"Found {article_count:,} articles on PubMed Central")

    # Combine texts
    if article_count == 0:
        combined_text = f'No articles found on PubMed Central for the query "{query}".'
    else:
        combined_text = format_references(
            references=tuple(texts),
            reference_type="paper",
            intro=f'Here are the top {article_count} articles on PubMed Central for the query "{query}":',
        )

    return combined_text


def run_retraction_check(ids: list[str]) -> str:
    """Checks papers for retraction notices and expressions of concern using the free SaneCite API.

    Only the IDs are sent. Notice data comes from the Retraction Watch database via Crossref.
    A failed lookup is reported as not checked, never as finding no notice.

    :param ids: Up to 10 DOIs, PubMed IDs, or PubMed Central IDs written with the PMC prefix.
    :return: A plain-text report for each paper.
    """
    print(f"Checking {len(ids)} papers for retraction notices")

    reports = []
    for index, paper_id in enumerate(ids[:10]):
        # The API allows 20 calls a minute
        if index > 0:
            time.sleep(3.1)
        try:
            response = requests.get("https://sanecite.com/api/check", params={"q": paper_id}, timeout=30)
            if response.status_code == 404:
                reports.append(f"{paper_id}: not found, so not checked.")
                continue
            response.raise_for_status()
            result = response.json()
            if result.get("ok") is not True or not isinstance(result.get("flagged"), list):
                raise ValueError("unexpected answer")
        except (requests.RequestException, ValueError):
            reports.append(f"{paper_id}: the check failed, so this paper was not checked.")
            continue

        lines = [f"{paper_id}:"]
        own = (result.get("paper") or {}).get("own")
        if own:
            lines.append(f"- This paper has a notice: {own.get('status')} ({own.get('notice_date')}).")
        else:
            lines.append("- No notice on record for this paper.")
        lines.append(f"- References checked: {result.get('refs_checked')} of {result.get('refs_total')}.")
        for notice in result["flagged"]:
            lines.append(
                f"- Cites {notice.get('doi') or notice.get('label')}: {notice.get('status')} ({notice.get('notice_date')})."
            )
        reports.append("\n".join(lines))

    if len(ids) > 10:
        reports.append(f"Only the first 10 of {len(ids)} papers were checked.")
    reports.append(
        "No notice on record does not mean a paper is sound. A cited work with a notice is not a finding against "
        "the citing paper, which may predate the notice or cite the work to criticize it. Source: Retraction Watch "
        "database, via Crossref, through the SaneCite API (https://sanecite.com/data)."
    )
    return "\n\n".join(reports)


def run_tools(
    tool_calls: list[ChatCompletionMessageToolCall],
) -> tuple[list[str], list[ChatCompletionMessageParam]]:
    """Runs the tools from chat completion tool calls.

    :param tool_calls: The tool calls from the chat completion response.
    :return: A tuple of (tool outputs as strings, tool response messages for API).
    """
    tool_outputs: list[str] = []
    tool_messages: list[ChatCompletionMessageParam] = []

    for tool_call in tool_calls:
        if tool_call.function.name == PUBMED_TOOL_NAME:
            # Extract the arguments from the tool call
            args_dict = json.loads(tool_call.function.arguments)

            # Run the tool
            output = run_pubmed_search(**args_dict)
            tool_outputs.append(output)

            # Create tool response message for the API
            tool_messages.append(
                {
                    "role": "tool",
                    "tool_call_id": tool_call.id,
                    "content": output,
                }
            )
        elif tool_call.function.name == RETRACTION_TOOL_NAME:
            args_dict = json.loads(tool_call.function.arguments)
            output = run_retraction_check(**args_dict)
            tool_outputs.append(output)
            tool_messages.append(
                {
                    "role": "tool",
                    "tool_call_id": tool_call.id,
                    "content": output,
                }
            )
        else:
            raise ValueError(f"Unknown tool: {tool_call.function.name}")

    return tool_outputs, tool_messages


def count_tokens(string: str, encoding_name: str = "cl100k_base") -> int:
    """Returns the number of tokens in a text string.

    :param string: The text string to count tokens in.
    :param encoding_name: The name of the encoding to use.
    :return: The number of tokens in the text string.
    """
    encoding = tiktoken.get_encoding(encoding_name)
    num_tokens = len(encoding.encode(string))

    return num_tokens


def update_token_counts(
    token_counts: dict[str, int],
    discussion: list[dict[str, str]],
    response: str,
) -> None:
    """Updates the token counts (in place) with a discussion and response.

    :param token_counts: The token counts to update.
    :param discussion: The discussion to update the token counts with.
    :param response: The response to update the token counts with.
    """
    new_input_token_count = sum(count_tokens(turn["message"]) for turn in discussion)
    new_output_token_count = count_tokens(response)

    token_counts["input"] += new_input_token_count
    token_counts["output"] += new_output_token_count

    token_counts["max"] = max(token_counts["max"], new_input_token_count + new_output_token_count)


def count_discussion_tokens(
    discussion: list[dict[str, str]],
) -> dict[str, int]:
    """Counts the number of tokens in a discussion.

    :param discussion: The discussion to count tokens in.
    :return: A dictionary of token counts.
    """
    token_counts = {
        "input": 0,
        "output": 0,
        "max": 0,
    }

    for index, turn in enumerate(discussion):
        if turn["agent"] != "User":
            update_token_counts(
                token_counts=token_counts,
                discussion=discussion[:index],
                response=turn["message"],
            )

    return token_counts


def _find_model_price_key(model: str, price_dict: dict[str, float]) -> str | None:
    """Finds the matching key in a price dictionary for a model.

    First checks for an exact match, then finds the longest prefix match.

    :param model: The name of the model.
    :param price_dict: The price dictionary to search.
    :return: The matching key or None if no match found.
    """
    if model in price_dict:
        return model

    # Find the longest prefix match
    matching_keys = [key for key in price_dict if model.startswith(key)]
    if matching_keys:
        return max(matching_keys, key=len)

    return None


def compute_token_cost(model: str, input_token_count: int, output_token_count: int) -> float:
    """Computes the token cost of a model given input and output token counts.

    :param model: The name of the model.
    :param input_token_count: The number of tokens in the input.
    :param output_token_count: The number of tokens in the output.
    :return: The token cost of the model.
    """
    input_key = _find_model_price_key(model, MODEL_TO_INPUT_PRICE_PER_TOKEN)
    output_key = _find_model_price_key(model, MODEL_TO_OUTPUT_PRICE_PER_TOKEN)

    if input_key is None or output_key is None:
        raise ValueError(f'Cost of model "{model}" not known')

    return (
        input_token_count * MODEL_TO_INPUT_PRICE_PER_TOKEN[input_key]
        + output_token_count * MODEL_TO_OUTPUT_PRICE_PER_TOKEN[output_key]
    )


def print_cost_and_time(
    token_counts: dict[str, int],
    model: str,
    elapsed_time: float,
) -> None:
    # Print token counts
    print(f"Input token count: {token_counts['input']:,}")
    print(f"Output token count: {token_counts['output']:,}")
    print(f"Tool token count: {token_counts['tool']:,}")
    print(f"Max token length: {token_counts['max']:,}")

    # Compute and print cost
    try:
        cost = compute_token_cost(
            model=model,
            input_token_count=token_counts["input"] + token_counts["tool"],
            output_token_count=token_counts["output"],
        )
        print(f"Cost: ${cost:.2f}")
    except ValueError as e:
        print(f"Warning: {e}")

    # Print time
    print(f"Time: {int(elapsed_time // 60)}:{int(elapsed_time % 60):02d}")


def compute_finetuning_cost(model: str, token_count: int, num_epochs: int = DEFAULT_FINETUNING_EPOCHS) -> float:
    """Computes the cost of fine-tuning a model.

    :param model: The model that will be finetuned.
    :param token_count: The number of training tokens for finetuning.
    :param num_epochs: Number of finetuning epochs.
    :return: The cost of finetuning.
    """
    if model not in FINETUNING_MODEL_TO_TRAINING_PRICE_PER_TOKEN:
        raise ValueError(f'Cost of model "{model}" not known')

    return token_count * FINETUNING_MODEL_TO_TRAINING_PRICE_PER_TOKEN[model] * num_epochs


def get_summary(discussion: list[dict[str, str]]) -> str:
    """Get the summary from a discussion.

    :param discussion: The discussion to extract the summary from.
    :return: The summary.
    """
    return discussion[-1]["message"]


def load_summaries(discussion_paths: list[Path]) -> tuple[str, ...]:
    """Load summaries from a list of discussion paths.

    :param discussion_paths: The paths to the discussion JSON files. The summary is the last entry in the discussion.
    :return: A tuple of summaries.
    """
    summaries = []
    for discussion_path in discussion_paths:
        with open(discussion_path, "r") as file:
            discussion = json.load(file)
        summaries.append(get_summary(discussion))

    return tuple(summaries)


def save_meeting(save_dir: Path, save_name: str, discussion: list[dict[str, str]]) -> None:
    """Save a meeting discussion to JSON and Markdown files.

    :param save_dir: The directory to save the discussion.
    :param save_name: The name of the discussion file that will be saved.
    :param discussion: The discussion to save.
    """
    # Create the save directory if it does not exist
    save_dir.mkdir(parents=True, exist_ok=True)

    # Save the discussion as JSON
    with open(save_dir / f"{save_name}.json", "w") as f:
        json.dump(discussion, f, indent=4)

    # Save the discussion as Markdown
    with open(save_dir / f"{save_name}.md", "w", encoding="utf-8") as file:
        for turn in discussion:
            file.write(f"## {turn['agent']}\n\n{turn['message']}\n\n")

"""Contains useful utility functions."""

import json
import urllib.parse
from pathlib import Path

import requests
import tiktoken
from openai import OpenAI
from openai.types.beta.threads.run import Run

from virtual_lab.constants import (
    MODEL_TO_INPUT_PRICE_PER_TOKEN,
    MODEL_TO_OUTPUT_PRICE_PER_TOKEN,
    PUBMED_TOOL_NAME,
)
from virtual_lab.prompts import format_references


def run_pubmed_search(
    query: str, num_articles: int = 3, abstract_only: bool = False, verbose: bool = True
) -> str:
    """Runs a PubMed search, returning the full text of the top matching article.

    :param query: The query to search PubMed with.
    :param num_articles: The number of articles to search for.
    :param abstract_only: Whether to return only the abstract.
    :param verbose: Whether to print verbose output.
    :return: The full text of the top matching article.
    """
    # Print search query
    if verbose:
        print(
            f'Searching PubMed Central for {num_articles} articles ({'abstracts' if abstract_only else 'full text'}) with query: "{query}"'
        )

    # Perform PubMed Central search for query to get PMC ID
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term={urllib.parse.quote_plus(query)}&retmax={2 * num_articles}&retmode=json&sort=relevance"
    response = requests.get(search_url)
    response.raise_for_status()
    pmc_ids_found = response.json()["esearchresult"]["idlist"]

    # Loop through top articles
    texts = []
    titles = []
    pmc_ids = []

    for pmc_id in pmc_ids_found:
        # Break if reached desired number of articles
        if len(pmc_ids) >= num_articles:
            break

        # Get full text of article from PMC ID in JSON form
        text_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_JSON/PMC{pmc_id}/unicode"
        response = requests.get(text_url)
        response.raise_for_status()

        try:
            article = response.json()
        except json.JSONDecodeError:
            if verbose:
                print(
                    f'Warning: Error decoding JSON for PMC ID {pmc_id}, "{response.text}"'
                )
            continue

        # Check if article is valid
        if len(article) != 1 or len(article[0]["documents"]) != 1:
            raise ValueError("Expected exactly one article")

        # Get document
        document = article[0]["documents"][0]

        # Get abstract or full text of article (excluding references)
        if abstract_only:
            passages = [
                passage
                for passage in document["passages"]
                if passage["infons"]["section_type"] in ["TITLE", "ABSTRACT"]
            ]
        else:
            passages = [
                passage
                for passage in document["passages"]
                if passage["infons"]["section_type"] != "REF"
            ]

        # Get text
        text = f"PMCID = {pmc_id}\n" + "\n".join(
            f"{passage['infons']['section_type']} ({passage['infons']['type']}): {passage['text']}"
            for passage in passages
        )

        # Get title
        title = next(
            passage["text"]
            for passage in passages
            if passage["infons"]["section_type"] == "TITLE"
        )

        texts.append(text)
        titles.append(title)
        pmc_ids.append(pmc_id)

    # Print articles found
    article_count = len(texts)

    if verbose:
        print(
            f"Found {article_count} articles on PubMed Central:\n"
            + "\n".join(f"{pmc_id}: {title}" for pmc_id, title in zip(pmc_ids, titles))
        )

    # Combine texts
    combined_text = format_references(
        references=tuple(texts),
        reference_type="paper",
        intro=f'Here are the top {article_count} articles on PubMed Central for the query "{query}":',
    )

    return combined_text


def run_tools(run: Run) -> list[dict[str, str]]:
    """Runs the tools in a required action.

    :param run: The run to run tools for.
    :return: A list of tool outputs.
    """
    # Define the list to store tool outputs
    tool_outputs = []

    # Loop through each tool in the required action and run it
    for tool in run.required_action.submit_tool_outputs.tool_calls:
        if tool.function.name == PUBMED_TOOL_NAME:
            # Extract the query from the tool arguments
            args_dict = json.loads(tool.function.arguments)

            # Run the tool and append the output to the list of tool outputs
            tool_outputs.append(
                {
                    "tool_call_id": tool.id,
                    "output": run_pubmed_search(**args_dict),
                }
            )
        else:
            raise ValueError(f"Unknown tool: {tool.function.name}")

    return tool_outputs


def get_messages(client: OpenAI, thread_id: str) -> list[dict]:
    """Gets messages from a thread.

    :param client: The OpenAI client.
    :param thread_id: The ID of the thread to get messages from.
    :return: A list of messages.
    """
    # Set up
    messages = []
    last_message = None
    params = {
        "thread_id": thread_id,
        "limit": 100,
        "order": "asc",
    }

    # Get all messages from the thread page by page
    while True:
        # Set up params
        if last_message is not None:
            params["after"] = last_message["id"]
        elif "after" in params:
            del params["after"]

        # Get page of messages
        new_messages = client.beta.threads.messages.list(**params)

        # Convert messages to dictionary
        new_messages: list[dict] = new_messages.to_dict()["data"]

        # Append new messages
        messages += new_messages

        # Break if no more messages
        if len(new_messages) < params["limit"]:
            break

        # Get last message
        last_message = messages[-1]

    # Verify all message content is length 1
    assert all(len(message["content"]) == 1 for message in messages)

    return messages


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

    token_counts["max"] = max(
        token_counts["max"], new_input_token_count + new_output_token_count
    )


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


def compute_token_cost(
    model: str, input_token_count: int, output_token_count: int
) -> float:
    """Computes the token cost of a model given input and output token counts.

    :param model: The name of the model.
    :param input_token_count: The number of tokens in the input.
    :param output_token_count: The number of tokens in the output.
    :return: The token cost of the model.
    """
    if (
        model not in MODEL_TO_INPUT_PRICE_PER_TOKEN
        or model not in MODEL_TO_OUTPUT_PRICE_PER_TOKEN
    ):
        raise ValueError(f'Cost of model "{model}" not known')

    return (
        input_token_count * MODEL_TO_INPUT_PRICE_PER_TOKEN[model]
        + output_token_count * MODEL_TO_OUTPUT_PRICE_PER_TOKEN[model]
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


def convert_messages_to_discussion(
    messages: list[dict], assistant_id_to_title: dict[str, str]
) -> list[dict[str, str]]:
    """Converts OpenAI messages into discussion format (list of message dictionaries).

    :param messages: The messages to convert.
    :param assistant_id_to_title: A dictionary mapping assistant IDs to titles.
    :return: The discussion format (list of message dictionaries).
    """
    return [
        {
            "agent": (
                assistant_id_to_title[message["assistant_id"]]
                if message["assistant_id"] is not None
                else "User"
            ),
            "message": message["content"][0]["text"]["value"],
        }
        for message in messages
    ]


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


def save_meeting(
    save_dir: Path, save_name: str, discussion: list[dict[str, str]]
) -> None:
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
    with open(save_dir / f"{save_name}.md", "w") as file:
        for turn in discussion:
            file.write(f"## {turn['agent']}\n\n{turn['message']}\n\n")

"""
Metabolic Target Finder using Genome-Scale Metabolic Models

This script implements comprehensive analysis tools for identifying gene
knockout/knockdown targets using constraint-based metabolic modeling.

Methods:
1. Single gene knockout analysis
2. Double gene knockout (synthetic lethality)
3. Essential gene identification
4. Growth-coupled production analysis
5. Flux variability analysis
6. Minimal media prediction

Author: Virtual Lab
Date: 2025-11-04
"""

import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import (
    flux_variability_analysis,
    single_gene_deletion,
    double_gene_deletion,
)
from cobra.sampling import sample
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Set
import warnings
warnings.filterwarnings('ignore')

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Argument parsing
from tap import Tap


class Args(Tap):
    """Argument parser for metabolic target finder"""

    model_file: str = None
    """Path to metabolic model file (SBML, JSON, or MAT format)"""

    model_id: str = "iML1515"
    """BiGG model ID to download if no file provided (e.g., iML1515, iJO1366, Recon3D)"""

    output_dir: str = "metabolic_results"
    """Directory to save results"""

    objective: str = None
    """Objective function (default: biomass reaction)"""

    target_metabolite: str = None
    """Target metabolite ID for production analysis"""

    ko_methods: List[str] = ["single", "essential", "fva"]
    """Methods to run: single, double, essential, synthetic_lethality, fva, sampling"""

    growth_threshold: float = 0.1
    """Minimum growth rate threshold (fraction of wildtype)"""

    n_samples: int = 1000
    """Number of flux samples for sampling analysis"""

    fva_fraction: float = 0.95
    """Fraction of optimal growth for FVA"""

    top_n_targets: int = 50
    """Number of top targets to report"""

    visualization: bool = True
    """Generate visualization plots"""

    random_seed: int = 42
    """Random seed for reproducibility"""


class MetabolicTargetFinder:
    """
    A comprehensive class for metabolic engineering target identification.

    This class uses constraint-based modeling with COBRApy to identify
    gene knockout/knockdown targets for various applications.
    """

    def __init__(self, args: Args):
        self.args = args
        self.model = None
        self.wildtype_growth = None
        self.results = {}

        # Set random seed
        np.random.seed(args.random_seed)

    def load_model(self) -> Model:
        """Load metabolic model from file or BiGG database"""
        print(f"Loading metabolic model...")

        if self.args.model_file:
            # Load from file
            print(f"  Loading from file: {self.args.model_file}")
            if self.args.model_file.endswith('.xml') or self.args.model_file.endswith('.sbml'):
                self.model = cobra.io.read_sbml_model(self.args.model_file)
            elif self.args.model_file.endswith('.json'):
                self.model = cobra.io.load_json_model(self.args.model_file)
            elif self.args.model_file.endswith('.mat'):
                self.model = cobra.io.load_matlab_model(self.args.model_file)
            else:
                raise ValueError(f"Unsupported file format: {self.args.model_file}")
        else:
            # Load from BiGG database
            print(f"  Downloading from BiGG database: {self.args.model_id}")
            try:
                self.model = cobra.io.load_model(self.args.model_id)
            except Exception as e:
                print(f"  Error: {e}")
                print(f"  Trying alternative method...")
                # Try to load E. coli core model as fallback
                self.model = cobra.io.load_model("textbook")
                print(f"  Loaded textbook E. coli core model as fallback")

        print(f"\n  Model: {self.model.id}")
        print(f"  Reactions: {len(self.model.reactions)}")
        print(f"  Metabolites: {len(self.model.metabolites)}")
        print(f"  Genes: {len(self.model.genes)}")

        # Set objective if specified
        if self.args.objective:
            self.model.objective = self.args.objective
            print(f"  Objective set to: {self.args.objective}")
        else:
            print(f"  Objective: {self.model.objective.expression}")

        # Calculate wildtype growth
        solution = self.model.optimize()
        self.wildtype_growth = solution.objective_value
        print(f"  Wildtype growth rate: {self.wildtype_growth:.4f}")

        return self.model

    def method_single_gene_knockout(self) -> pd.DataFrame:
        """
        Perform single gene knockout analysis.

        Tests the effect of deleting each gene individually on growth rate.

        Returns:
            DataFrame with gene IDs, growth rates, and growth fractions
        """
        print("\n=== Single Gene Knockout Analysis ===")
        print(f"Analyzing {len(self.model.genes)} genes...")

        # Perform single gene deletion
        deletion_results = single_gene_deletion(self.model)

        # Process results
        results = []
        for idx, row in deletion_results.iterrows():
            gene_id = row['ids']
            growth = row['growth']
            growth_fraction = growth / self.wildtype_growth if self.wildtype_growth > 0 else 0

            results.append({
                'gene_id': gene_id,
                'growth_rate': growth,
                'growth_fraction': growth_fraction,
                'growth_reduction': 1 - growth_fraction,
                'status': row['status']
            })

        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('growth_reduction', ascending=False)

        # Categorize genes
        essential = results_df[results_df['growth_fraction'] < self.args.growth_threshold]
        nonessential = results_df[results_df['growth_fraction'] >= self.args.growth_threshold]

        print(f"\nResults:")
        print(f"  Essential genes: {len(essential)} ({len(essential)/len(results_df)*100:.1f}%)")
        print(f"  Non-essential genes: {len(nonessential)} ({len(nonessential)/len(results_df)*100:.1f}%)")
        print(f"  Growth-reducing genes: {(results_df['growth_reduction'] > 0.1).sum()}")

        self.results['single_knockout'] = results_df
        return results_df

    def method_essential_genes(self) -> pd.DataFrame:
        """
        Identify essential genes.

        Essential genes are those whose deletion results in growth below threshold.

        Returns:
            DataFrame with essential genes and their properties
        """
        print("\n=== Essential Gene Analysis ===")

        if 'single_knockout' not in self.results:
            print("Running single gene knockout first...")
            self.method_single_gene_knockout()

        results_df = self.results['single_knockout']
        essential = results_df[results_df['growth_fraction'] < self.args.growth_threshold].copy()

        # Add gene information
        essential_with_info = []
        for _, row in essential.iterrows():
            gene_id = row['gene_id']
            gene = self.model.genes.get_by_id(gene_id)

            # Get reactions involving this gene
            reactions = list(gene.reactions)
            reaction_ids = [r.id for r in reactions]

            essential_with_info.append({
                'gene_id': gene_id,
                'gene_name': gene.name,
                'growth_rate': row['growth_rate'],
                'growth_fraction': row['growth_fraction'],
                'n_reactions': len(reactions),
                'reaction_ids': '; '.join(reaction_ids[:5])  # First 5 reactions
            })

        essential_df = pd.DataFrame(essential_with_info)
        essential_df = essential_df.sort_values('growth_fraction')

        print(f"\nIdentified {len(essential_df)} essential genes")
        print(f"\nTop 10 most critical essential genes:")
        print(essential_df.head(10)[['gene_id', 'gene_name', 'growth_fraction', 'n_reactions']].to_string(index=False))

        self.results['essential_genes'] = essential_df
        return essential_df

    def method_double_gene_knockout(self, gene_list: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Perform double gene knockout analysis (synthetic lethality).

        Args:
            gene_list: List of gene IDs to test. If None, uses top non-essential genes.

        Returns:
            DataFrame with gene pair knockouts and growth effects
        """
        print("\n=== Double Gene Knockout Analysis (Synthetic Lethality) ===")

        if gene_list is None:
            # Use top non-essential genes that have some growth effect
            if 'single_knockout' not in self.results:
                self.method_single_gene_knockout()

            single_ko = self.results['single_knockout']
            # Select genes that are non-essential but reduce growth
            candidates = single_ko[
                (single_ko['growth_fraction'] >= self.args.growth_threshold) &
                (single_ko['growth_reduction'] > 0.05)
            ].head(50)  # Top 50 candidates

            gene_list = candidates['gene_id'].tolist()

        print(f"Testing {len(gene_list)} genes...")
        print(f"Total combinations: {len(gene_list) * (len(gene_list) - 1) // 2}")

        # Perform double gene deletion
        try:
            deletion_results = double_gene_deletion(
                self.model,
                gene_list1=gene_list,
                gene_list2=gene_list
            )

            # Process results
            results = []
            for idx, row in deletion_results.iterrows():
                gene_ids = row['ids'].split(', ')
                if len(gene_ids) == 2 and gene_ids[0] != gene_ids[1]:  # Exclude self-pairs
                    growth = row['growth']
                    growth_fraction = growth / self.wildtype_growth if self.wildtype_growth > 0 else 0

                    # Get single knockout growth for comparison
                    single_ko = self.results['single_knockout']
                    gene1_growth = single_ko[single_ko['gene_id'] == gene_ids[0]]['growth_fraction'].iloc[0]
                    gene2_growth = single_ko[single_ko['gene_id'] == gene_ids[1]]['growth_fraction'].iloc[0]

                    # Synergy score: negative means synthetic lethality
                    expected_growth = min(gene1_growth, gene2_growth)
                    synergy = growth_fraction - expected_growth

                    results.append({
                        'gene1': gene_ids[0],
                        'gene2': gene_ids[1],
                        'double_ko_growth': growth,
                        'double_ko_fraction': growth_fraction,
                        'gene1_single_fraction': gene1_growth,
                        'gene2_single_fraction': gene2_growth,
                        'expected_growth': expected_growth,
                        'synergy_score': synergy,
                        'synthetic_lethal': growth_fraction < self.args.growth_threshold and
                                          gene1_growth >= self.args.growth_threshold and
                                          gene2_growth >= self.args.growth_threshold,
                        'status': row['status']
                    })

            results_df = pd.DataFrame(results)
            results_df = results_df.sort_values('synergy_score')

            # Find synthetic lethal pairs
            synthetic_lethal = results_df[results_df['synthetic_lethal'] == True]

            print(f"\nResults:")
            print(f"  Total pairs tested: {len(results_df)}")
            print(f"  Synthetic lethal pairs: {len(synthetic_lethal)}")
            if len(synthetic_lethal) > 0:
                print(f"\nTop 10 synthetic lethal pairs:")
                print(synthetic_lethal.head(10)[['gene1', 'gene2', 'double_ko_fraction', 'synergy_score']].to_string(index=False))

            self.results['double_knockout'] = results_df
            return results_df

        except Exception as e:
            print(f"Error in double knockout analysis: {e}")
            print("Skipping double knockout analysis")
            return pd.DataFrame()

    def method_flux_variability_analysis(self) -> pd.DataFrame:
        """
        Perform flux variability analysis (FVA).

        FVA identifies the range of possible flux values for each reaction
        while maintaining a certain fraction of optimal growth.

        Returns:
            DataFrame with flux ranges for all reactions
        """
        print("\n=== Flux Variability Analysis (FVA) ===")
        print(f"FVA at {self.args.fva_fraction*100}% of optimal growth...")

        # Perform FVA
        fva_results = flux_variability_analysis(
            self.model,
            fraction_of_optimum=self.args.fva_fraction
        )

        # Calculate flux spans
        fva_results['flux_span'] = fva_results['maximum'] - fva_results['minimum']
        fva_results['flux_mid'] = (fva_results['maximum'] + fva_results['minimum']) / 2
        fva_results['flexibility'] = fva_results['flux_span'] / (abs(fva_results['flux_mid']) + 1e-9)

        # Add reaction information
        reaction_info = []
        for rxn_id in fva_results.index:
            rxn = self.model.reactions.get_by_id(rxn_id)
            genes = [g.id for g in rxn.genes]

            reaction_info.append({
                'reaction_id': rxn_id,
                'reaction_name': rxn.name,
                'n_genes': len(genes),
                'gene_ids': '; '.join(genes[:3]) if genes else 'No genes',
                'subsystem': rxn.subsystem
            })

        info_df = pd.DataFrame(reaction_info).set_index('reaction_id')
        fva_results = fva_results.join(info_df)

        # Identify reactions with zero flux (non-essential)
        zero_flux = fva_results[(fva_results['minimum'] == 0) & (fva_results['maximum'] == 0)]
        flexible = fva_results[fva_results['flux_span'] > 1.0]
        essential = fva_results[(fva_results['minimum'] > 0) | (fva_results['maximum'] < 0)]

        print(f"\nResults:")
        print(f"  Total reactions: {len(fva_results)}")
        print(f"  Zero flux (potentially non-essential): {len(zero_flux)}")
        print(f"  Highly flexible reactions: {len(flexible)}")
        print(f"  Essential flux reactions: {len(essential)}")

        self.results['fva'] = fva_results.sort_values('flux_span', ascending=False)
        return self.results['fva']

    def method_growth_coupled_production(self, target_metabolite: str) -> pd.DataFrame:
        """
        Identify gene knockout targets for growth-coupled production.

        Args:
            target_metabolite: Metabolite ID to produce

        Returns:
            DataFrame with gene knockouts and production-to-growth ratios
        """
        print(f"\n=== Growth-Coupled Production Analysis ===")
        print(f"Target metabolite: {target_metabolite}")

        if 'single_knockout' not in self.results:
            self.method_single_gene_knockout()

        # Get metabolite
        try:
            metabolite = self.model.metabolites.get_by_id(target_metabolite)
        except KeyError:
            print(f"Error: Metabolite {target_metabolite} not found in model")
            return pd.DataFrame()

        # Create exchange reaction if not exists
        exchange_rxn_id = f"EX_{target_metabolite}"
        if exchange_rxn_id not in self.model.reactions:
            print(f"Creating exchange reaction: {exchange_rxn_id}")
            exchange = Reaction(exchange_rxn_id)
            exchange.name = f"Exchange for {metabolite.name}"
            exchange.lower_bound = 0
            exchange.upper_bound = 1000
            exchange.add_metabolites({metabolite: -1})
            self.model.add_reactions([exchange])

        # Test each knockout for production coupling
        single_ko = self.results['single_knockout']
        viable_ko = single_ko[single_ko['growth_fraction'] >= self.args.growth_threshold]

        results = []
        for _, row in viable_ko.iterrows():
            gene_id = row['gene_id']
            gene = self.model.genes.get_by_id(gene_id)

            # Perform knockout
            with self.model:
                gene.knock_out()

                # Optimize for growth
                solution = self.model.optimize()
                ko_growth = solution.objective_value

                if ko_growth >= self.wildtype_growth * self.args.growth_threshold:
                    # Get production flux
                    production = solution.fluxes[exchange_rxn_id]

                    # Calculate production-to-growth ratio
                    ratio = production / ko_growth if ko_growth > 0 else 0

                    results.append({
                        'gene_id': gene_id,
                        'ko_growth': ko_growth,
                        'production': production,
                        'production_per_growth': ratio,
                        'growth_fraction': ko_growth / self.wildtype_growth
                    })

        results_df = pd.DataFrame(results)
        if len(results_df) > 0:
            results_df = results_df.sort_values('production_per_growth', ascending=False)

            print(f"\nResults:")
            print(f"  Viable knockouts tested: {len(viable_ko)}")
            print(f"  Knockouts with production: {(results_df['production'] > 0).sum()}")
            print(f"\nTop 10 production-coupled knockouts:")
            print(results_df.head(10).to_string(index=False))

            self.results['growth_coupled'] = results_df
        else:
            print("No viable production knockouts found")

        return results_df

    def method_flux_sampling(self) -> pd.DataFrame:
        """
        Perform flux sampling to explore solution space.

        Returns:
            DataFrame with flux statistics across samples
        """
        print(f"\n=== Flux Sampling Analysis ===")
        print(f"Generating {self.args.n_samples} flux samples...")

        try:
            # Sample flux space
            samples = sample(self.model, self.args.n_samples, seed=self.args.random_seed)

            # Calculate statistics
            flux_stats = pd.DataFrame({
                'mean_flux': samples.mean(),
                'std_flux': samples.std(),
                'min_flux': samples.min(),
                'max_flux': samples.max(),
                'cv': samples.std() / (samples.mean().abs() + 1e-9)  # Coefficient of variation
            })

            # Add reaction information
            reaction_info = []
            for rxn_id in flux_stats.index:
                rxn = self.model.reactions.get_by_id(rxn_id)
                genes = [g.id for g in rxn.genes]

                reaction_info.append({
                    'reaction_id': rxn_id,
                    'reaction_name': rxn.name,
                    'n_genes': len(genes),
                    'gene_ids': '; '.join(genes[:3]) if genes else 'No genes',
                    'subsystem': rxn.subsystem
                })

            info_df = pd.DataFrame(reaction_info).set_index('reaction_id')
            flux_stats = flux_stats.join(info_df)

            # Identify highly variable reactions
            high_cv = flux_stats[flux_stats['cv'] > 1.0].sort_values('cv', ascending=False)

            print(f"\nResults:")
            print(f"  Reactions sampled: {len(flux_stats)}")
            print(f"  High variability reactions (CV > 1): {len(high_cv)}")

            self.results['flux_sampling'] = flux_stats
            return flux_stats

        except Exception as e:
            print(f"Error in flux sampling: {e}")
            print("Flux sampling requires additional dependencies (optlang)")
            return pd.DataFrame()

    def identify_best_targets(self, criteria: str = "growth_reduction") -> pd.DataFrame:
        """
        Identify best knockout targets based on specified criteria.

        Args:
            criteria: Selection criteria (growth_reduction, essential, synthetic_lethal)

        Returns:
            DataFrame with ranked targets
        """
        print(f"\n=== Identifying Best Targets ({criteria}) ===")

        if criteria == "growth_reduction":
            if 'single_knockout' not in self.results:
                self.method_single_gene_knockout()

            results = self.results['single_knockout']
            # Non-essential genes with significant growth reduction
            targets = results[
                (results['growth_fraction'] >= self.args.growth_threshold) &
                (results['growth_reduction'] > 0.1)
            ]
            targets = targets.sort_values('growth_reduction', ascending=False)

        elif criteria == "essential":
            if 'essential_genes' not in self.results:
                self.method_essential_genes()
            targets = self.results['essential_genes']

        elif criteria == "synthetic_lethal":
            if 'double_knockout' not in self.results:
                self.method_double_gene_knockout()
            targets = self.results['double_knockout']
            targets = targets[targets['synthetic_lethal'] == True]

        else:
            raise ValueError(f"Unknown criteria: {criteria}")

        print(f"\nTop {min(self.args.top_n_targets, len(targets))} targets:")
        print(targets.head(self.args.top_n_targets).to_string(index=False))

        return targets

    def visualize_results(self, output_dir: str):
        """Generate visualization plots for the results"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        pdf_path = os.path.join(output_dir, 'metabolic_analysis.pdf')

        with PdfPages(pdf_path) as pdf:
            # Plot 1: Single knockout growth distribution
            if 'single_knockout' in self.results:
                self._plot_knockout_distribution(pdf)

            # Plot 2: Essential genes by category
            if 'essential_genes' in self.results:
                self._plot_essential_genes(pdf)

            # Plot 3: Synthetic lethality network
            if 'double_knockout' in self.results:
                self._plot_synthetic_lethality(pdf)

            # Plot 4: FVA flux ranges
            if 'fva' in self.results:
                self._plot_fva_results(pdf)

            # Plot 5: Growth-coupled production
            if 'growth_coupled' in self.results:
                self._plot_production_coupling(pdf)

        print(f"\nVisualizations saved to {pdf_path}")

    def _plot_knockout_distribution(self, pdf):
        """Plot distribution of growth rates after single knockout"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        df = self.results['single_knockout']

        # Histogram
        ax = axes[0]
        ax.hist(df['growth_fraction'], bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(self.args.growth_threshold, color='red', linestyle='--',
                  label=f'Essential threshold ({self.args.growth_threshold})')
        ax.set_xlabel('Growth Fraction (vs wildtype)', fontsize=12)
        ax.set_ylabel('Number of Genes', fontsize=12)
        ax.set_title('Distribution of Single Gene Knockout Effects', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)

        # Pie chart
        ax = axes[1]
        essential = (df['growth_fraction'] < self.args.growth_threshold).sum()
        reducing = ((df['growth_fraction'] >= self.args.growth_threshold) &
                   (df['growth_reduction'] > 0.1)).sum()
        minimal = len(df) - essential - reducing

        sizes = [essential, reducing, minimal]
        labels = [f'Essential\n({essential})',
                 f'Growth-reducing\n({reducing})',
                 f'Minimal effect\n({minimal})']
        colors = ['#ff6b6b', '#ffd93d', '#6bcf7f']

        ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
              startangle=90, textprops={'fontsize': 11})
        ax.set_title('Gene Classification by Knockout Effect', fontsize=14, fontweight='bold')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    def _plot_essential_genes(self, pdf):
        """Plot essential genes analysis"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        df = self.results['essential_genes'].head(30)

        # Bar plot
        ax = axes[0]
        bars = ax.barh(range(len(df)), df['growth_fraction'])
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['gene_id'], fontsize=8)
        ax.set_xlabel('Growth Fraction', fontsize=12)
        ax.set_title('Top 30 Essential Genes', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Reactions per gene
        ax = axes[1]
        ax.scatter(df['n_reactions'], df['growth_fraction'], s=100, alpha=0.6)
        ax.set_xlabel('Number of Reactions', fontsize=12)
        ax.set_ylabel('Growth Fraction', fontsize=12)
        ax.set_title('Essential Genes: Reactions vs Growth', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    def _plot_synthetic_lethality(self, pdf):
        """Plot synthetic lethality results"""
        if 'double_knockout' not in self.results or len(self.results['double_knockout']) == 0:
            return

        df = self.results['double_knockout']

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Scatter plot of synergy scores
        ax = axes[0]
        synthetic = df[df['synthetic_lethal'] == True]
        non_synthetic = df[df['synthetic_lethal'] == False]

        ax.scatter(non_synthetic['expected_growth'], non_synthetic['double_ko_fraction'],
                  alpha=0.3, s=20, color='gray', label='Non-synthetic')
        ax.scatter(synthetic['expected_growth'], synthetic['double_ko_fraction'],
                  alpha=0.7, s=50, color='red', label='Synthetic lethal')
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Expected')
        ax.axhline(self.args.growth_threshold, color='blue', linestyle='--', alpha=0.5)
        ax.set_xlabel('Expected Growth (min of single KOs)', fontsize=12)
        ax.set_ylabel('Observed Double KO Growth', fontsize=12)
        ax.set_title('Synthetic Lethality Analysis', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)

        # Histogram of synergy scores
        ax = axes[1]
        ax.hist(df['synergy_score'], bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(0, color='red', linestyle='--', label='No synergy')
        ax.set_xlabel('Synergy Score', fontsize=12)
        ax.set_ylabel('Number of Gene Pairs', fontsize=12)
        ax.set_title('Distribution of Synergy Scores', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    def _plot_fva_results(self, pdf):
        """Plot FVA flux variability results"""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        df = self.results['fva']

        # Flux span distribution
        ax = axes[0, 0]
        ax.hist(df['flux_span'], bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Flux Span (max - min)', fontsize=12)
        ax.set_ylabel('Number of Reactions', fontsize=12)
        ax.set_title('FVA: Flux Span Distribution', fontsize=13, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(alpha=0.3)

        # Min vs Max flux
        ax = axes[0, 1]
        ax.scatter(df['minimum'], df['maximum'], alpha=0.5, s=10)
        ax.plot([-100, 100], [-100, 100], 'r--', alpha=0.5)
        ax.set_xlabel('Minimum Flux', fontsize=12)
        ax.set_ylabel('Maximum Flux', fontsize=12)
        ax.set_title('FVA: Min vs Max Flux', fontsize=13, fontweight='bold')
        ax.grid(alpha=0.3)

        # Top flexible reactions
        ax = axes[1, 0]
        top_flex = df.nlargest(20, 'flux_span')
        bars = ax.barh(range(len(top_flex)), top_flex['flux_span'])
        ax.set_yticks(range(len(top_flex)))
        ax.set_yticklabels(top_flex.index, fontsize=8)
        ax.set_xlabel('Flux Span', fontsize=12)
        ax.set_title('Top 20 Most Flexible Reactions', fontsize=13, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Subsystem analysis
        ax = axes[1, 1]
        subsys_flex = df.groupby('subsystem')['flux_span'].mean().sort_values(ascending=False).head(15)
        bars = ax.barh(range(len(subsys_flex)), subsys_flex.values)
        ax.set_yticks(range(len(subsys_flex)))
        ax.set_yticklabels(subsys_flex.index, fontsize=9)
        ax.set_xlabel('Average Flux Span', fontsize=12)
        ax.set_title('Average Flexibility by Subsystem', fontsize=13, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    def _plot_production_coupling(self, pdf):
        """Plot growth-coupled production results"""
        df = self.results['growth_coupled']

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Production vs Growth
        ax = axes[0]
        ax.scatter(df['ko_growth'], df['production'], s=100, alpha=0.6)
        ax.set_xlabel('Growth Rate (KO)', fontsize=12)
        ax.set_ylabel('Production Flux', fontsize=12)
        ax.set_title('Growth-Coupled Production', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3)

        # Top targets
        ax = axes[1]
        top = df.nlargest(20, 'production_per_growth')
        bars = ax.barh(range(len(top)), top['production_per_growth'])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top['gene_id'], fontsize=9)
        ax.set_xlabel('Production per Growth', fontsize=12)
        ax.set_title('Top 20 Production-Coupled Knockouts', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    def save_results(self, output_dir: str):
        """Save all results to CSV files"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        # Save individual method results
        for method_name, results_df in self.results.items():
            if len(results_df) > 0:
                output_path = os.path.join(output_dir, f'{method_name}_results.csv')
                results_df.to_csv(output_path, index=False if isinstance(results_df.index, pd.RangeIndex) else True)
                print(f"Saved {method_name} results to {output_path}")

        # Save model summary
        summary = {
            'model_id': self.model.id,
            'n_reactions': len(self.model.reactions),
            'n_metabolites': len(self.model.metabolites),
            'n_genes': len(self.model.genes),
            'wildtype_growth': self.wildtype_growth,
            'objective': str(self.model.objective.expression)
        }

        summary_df = pd.DataFrame([summary])
        summary_path = os.path.join(output_dir, 'model_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"Saved model summary to {summary_path}")


def main():
    """Main execution function"""
    # Parse arguments
    args = Args().parse_args()

    print("="*80)
    print("METABOLIC TARGET FINDER - GENOME-SCALE METABOLIC MODELING")
    print("="*80)
    print(f"\nConfiguration:")
    if args.model_file:
        print(f"  Model file: {args.model_file}")
    else:
        print(f"  Model ID: {args.model_id}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Methods: {', '.join(args.ko_methods)}")
    print(f"  Growth threshold: {args.growth_threshold}")
    print()

    # Initialize finder
    finder = MetabolicTargetFinder(args)

    # Load model
    finder.load_model()

    # Run selected methods
    if 'single' in args.ko_methods:
        finder.method_single_gene_knockout()

    if 'essential' in args.ko_methods:
        finder.method_essential_genes()

    if 'double' in args.ko_methods:
        finder.method_double_gene_knockout()

    if 'fva' in args.ko_methods:
        finder.method_flux_variability_analysis()

    if 'sampling' in args.ko_methods:
        finder.method_flux_sampling()

    if 'production' in args.ko_methods and args.target_metabolite:
        finder.method_growth_coupled_production(args.target_metabolite)

    # Identify best targets
    if 'single_knockout' in finder.results:
        finder.identify_best_targets(criteria="growth_reduction")

    # Save results
    finder.save_results(args.output_dir)

    # Generate visualizations
    if args.visualization and len(finder.results) > 0:
        print("\n=== Generating Visualizations ===")
        finder.visualize_results(args.output_dir)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults saved to: {args.output_dir}/")
    print(f"  - Model summary: model_summary.csv")
    for method in finder.results.keys():
        print(f"  - {method}: {method}_results.csv")
    if args.visualization:
        print(f"  - Visualizations: metabolic_analysis.pdf")
    print()


if __name__ == "__main__":
    main()

"""
Example: Simple Metabolic Modeling Analysis

This script demonstrates a simple metabolic modeling workflow
using the MetabolicTargetFinder class.
"""

from metabolic_target_finder import MetabolicTargetFinder, Args

def main():
    """Run a simple metabolic modeling example"""

    print("="*80)
    print("SIMPLE METABOLIC MODELING EXAMPLE")
    print("="*80)
    print()

    # Configure analysis
    args = Args(
        model_id="textbook",  # E. coli core model (fast)
        output_dir="../../example_metabolic_results",
        ko_methods=["single", "essential", "fva"],
        growth_threshold=0.1,
        visualization=True,
        random_seed=42
    )

    print("Configuration:")
    print(f"  Model: {args.model_id}")
    print(f"  Output: {args.output_dir}")
    print(f"  Methods: {', '.join(args.ko_methods)}")
    print(f"  Growth threshold: {args.growth_threshold}")
    print()

    # Initialize finder
    finder = MetabolicTargetFinder(args)

    # Load model
    print("Step 1: Loading metabolic model...")
    model = finder.load_model()
    print(f"  ✓ Loaded {model.id} with {len(model.genes)} genes")
    print(f"  ✓ Wildtype growth: {finder.wildtype_growth:.4f}")
    print()

    # Run single gene knockout
    print("Step 2: Single gene knockout analysis...")
    single_ko = finder.method_single_gene_knockout()
    print(f"  ✓ Analyzed {len(single_ko)} genes")
    print()

    # Identify essential genes
    print("Step 3: Essential gene identification...")
    essential = finder.method_essential_genes()
    print(f"  ✓ Found {len(essential)} essential genes")
    print()

    # Flux variability analysis
    print("Step 4: Flux variability analysis...")
    fva = finder.method_flux_variability_analysis()
    print(f"  ✓ Analyzed {len(fva)} reactions")
    print()

    # Identify best targets
    print("Step 5: Identifying best targets...")
    targets = finder.identify_best_targets(criteria="growth_reduction")
    print(f"  ✓ Found {len(targets)} potential targets")
    print()

    # Display top targets
    print("="*80)
    print("TOP 10 GENE KNOCKOUT TARGETS")
    print("="*80)
    print()
    print(targets.head(10)[['gene_id', 'growth_rate', 'growth_fraction', 'growth_reduction']].to_string(index=False))
    print()

    # Save results
    print("="*80)
    print("Step 6: Saving results...")
    finder.save_results(args.output_dir)
    print(f"  ✓ Results saved to {args.output_dir}/")
    print()

    # Generate visualizations
    print("Step 7: Generating visualizations...")
    finder.visualize_results(args.output_dir)
    print(f"  ✓ Visualizations saved to {args.output_dir}/metabolic_analysis.pdf")
    print()

    print("="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print()
    print("Next steps:")
    print("  1. Review results in:", args.output_dir)
    print("  2. Examine visualizations in: metabolic_analysis.pdf")
    print("  3. Try with larger models (iML1515, Recon3D)")
    print("  4. Explore double knockout and production analysis")
    print()


if __name__ == "__main__":
    main()

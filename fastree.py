import tempfile
from contextlib import contextmanager

FASTREER_VERSION = "2.1.3"
FASTREER_VERSION = "2.1.4"

# Default paths for BioFM resources (can be overridden via environment variables or CLI args)
DEFAULT_BIOFM_MODEL = os.environ.get("BIOFM_MODEL", "m42-health/BioFM-265M")
@@ -416,17 +416,35 @@ def add_embedding_args(p):
                       choices=["CHROM_POS", "CHROM_POS_REF_ALT", "VCF_ID"],
                       help="Variant key format for embedding lookup (default: CHROM_POS_REF_ALT)")

    def add_window_args(p):
        p.add_argument("--window-bp", dest="window_bp", type=int, default=None,
                       help="Emit one matrix/tree per window of N base pairs "
                            "(mutually exclusive with --window-variants)")
        p.add_argument("--window-variants", dest="window_variants", type=int, default=None,
                       help="Emit one matrix/tree per N consecutive variants "
                            "(mutually exclusive with --window-bp)")
        p.add_argument("--step", dest="window_step", type=int, default=None,
                       help="Window step (defaults to window size, i.e. tiled). "
                            "Sliding windows not yet implemented.")
        p.add_argument("--min-variants", dest="window_min_variants", type=int, default=None,
                       help="Minimum number of variants required to emit a window (default 1)")
    # Subcommand for VCF-based distance matrix
    parser_vcf2dist = subparsers.add_parser("VCF2DIST", help="Compute distance matrix from VCF(s)")
    add_common_input_output(parser_vcf2dist)
    add_common_vcf_args(parser_vcf2dist)
    add_embedding_args(parser_vcf2dist)
    add_window_args(parser_vcf2dist)
    parser_vcf2dist.add_argument("--long", dest="long_format", action="store_true",
                                 help="Emit long-form TSV (chrom, start, end, sample_i, sample_j, dist) "
                                      "instead of concatenated matrices (windowed mode only)")

    # Subcommand for VCF-based tree
    parser_vcf2tree = subparsers.add_parser("VCF2TREE", help="Compute tree from VCF(s)")
    add_common_input_output(parser_vcf2tree)
    add_common_vcf_args(parser_vcf2tree)
    add_embedding_args(parser_vcf2tree)
    add_window_args(parser_vcf2tree)
    parser_vcf2tree.add_argument("-b", "--bootstrap", type=int, default=0,
                                 help="Number of bootstrap replicates to perform (default: 0, no bootstrapping)")

@@ -540,6 +558,22 @@ def resolve_inputs(args, allow_multiple=True):
            params.extend(["--embeddings-format", args.embeddings_format])
        if getattr(args, 'variant_key', None) and args.variant_key != "CHROM_POS_REF_ALT":
            params.extend(["--variant-key", args.variant_key])
        # forward windowing options if provided
        window_bp = getattr(args, 'window_bp', None)
        window_variants = getattr(args, 'window_variants', None)
        if window_bp is not None and window_variants is not None:
            print("Error: --window-bp and --window-variants are mutually exclusive.", file=sys.stderr)
            sys.exit(1)
        if window_bp is not None:
            params.extend(["--window-bp", str(int(window_bp))])
        if window_variants is not None:
            params.extend(["--window-variants", str(int(window_variants))])
        if getattr(args, 'window_step', None) is not None:
            params.extend(["--step", str(int(args.window_step))])
        if getattr(args, 'window_min_variants', None) is not None:
            params.extend(["--min-variants", str(int(args.window_min_variants))])
        if args.command == "VCF2DIST" and getattr(args, 'long_format', False):
            params.append("--long")
        for f in input_files:
            params.extend(["-i", f])
        run_java_tool(args.command, params, args.lib, args.mem, args.output, args.verbose, args.extraVerbose, args.pipe_stderr)

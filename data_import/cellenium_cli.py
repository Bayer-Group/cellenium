import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", dest="subcommand")
    parser_study_import = subparsers.add_parser("study-import", help="cellenium study import tool")
    parser_study_import.add_argument(
        "filename",
        help="h5ad/h5mu file created for cellenium (e.g. using a jupyter lab notebook), local or S3 file",
        type=str,
    )
    parser_study_import.add_argument(
        "--analyze-database",
        help="analyses the database schema after insert of study",
        action="store_true",
    )

    parser_deg_calc = subparsers.add_parser(
        "differential-expression-calculation", help="calculate differential expression for a study that is already imported in cellenium"
    )
    parser_deg_calc.add_argument(
        "study_id",
        help="study ID for which differential expressed genes are calculated",
        type=int,
    )
    parser_deg_calc.add_argument(
        "annotation_group_ids",
        help="annotation group ID to recalculate, e.g. the ID for celltype; can be multiple IDs",
        type=int,
        nargs="+",
    )
    args = parser.parse_args()

    if args.subcommand == "study-import":
        import study_import

        study_import.import_study(args.filename, args.analyze_database)
    if args.subcommand == "differential-expression-calculation":
        import differential_expression_calculation

        differential_expression_calculation.differential_expression_calculation(args.study_id, args.annotation_group_ids)

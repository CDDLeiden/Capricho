"""Explore the ChEMBL database schema without loading the entire database."""

import pandas as pd
from tabulate import tabulate

import chembl_downloader

from .downloader import check_and_download_chembl_db


def get_tables(conn):
    """Get a list of all tables in the database."""
    query = """
    SELECT name FROM sqlite_master
    WHERE type='table'
    ORDER BY name;
    """
    return pd.read_sql(query, conn)


def get_table_info(conn, table_name):
    """Get detailed information about a specific table."""
    query = f"PRAGMA table_info({table_name});"
    return pd.read_sql(query, conn)


def get_table_counts(conn, table_name):
    """Get the row count for a specific table."""
    query = f"SELECT COUNT(*) as count FROM {table_name};"
    return pd.read_sql(query, conn)["count"][0]


def get_sample_data(conn, table_name, limit=5):
    """Get sample data from a table."""
    query = f"SELECT * FROM {table_name} LIMIT {limit};"
    return pd.read_sql(query, conn)


def get_foreign_keys(conn, table_name):
    """Get foreign key information for a table."""
    query = f"PRAGMA foreign_key_list({table_name});"
    return pd.read_sql(query, conn)


def find_related_tables(conn, table_name):
    """Find tables that reference this table or are referenced by this table."""
    # Get all tables
    tables = get_tables(conn)["name"].tolist()

    # Find relationships
    relationships = []

    # Look for foreign keys in the current table
    fks_in_table = get_foreign_keys(conn, table_name)
    for _, row in fks_in_table.iterrows():
        relationships.append(
            {
                "type": "outgoing",
                "table": table_name,
                "column": row["from"],
                "references_table": row["table"],
                "references_column": row["to"],
            }
        )

    # Look for tables that reference this table
    for other_table in tables:
        if other_table == table_name:
            continue

        fks = get_foreign_keys(conn, other_table)
        for _, row in fks.iterrows():
            if row["table"] == table_name:
                relationships.append(
                    {
                        "type": "incoming",
                        "table": other_table,
                        "column": row["from"],
                        "references_table": table_name,
                        "references_column": row["to"],
                    }
                )

    return relationships


def search_tables_for_column(conn, column_pattern):
    """Search for tables that have columns matching the given pattern."""
    tables = get_tables(conn)["name"].tolist()
    results = []

    for table in tables:
        columns = get_table_info(conn, table)
        matches = columns[columns["name"].str.contains(column_pattern, case=False)]

        if not matches.empty:
            for _, match in matches.iterrows():
                results.append({"table": table, "column": match["name"], "type": match["type"]})

    return results


def list_all_tables(conn):
    """List all tables with their row counts."""
    tables = get_tables(conn)

    results = []
    for table_name in tables["name"]:
        try:
            count = get_table_counts(conn, table_name)
            results.append({"table": table_name, "rows": count})
        except Exception as e:
            results.append({"table": table_name, "rows": f"Error: {str(e)}"})

    return pd.DataFrame(results).sort_values(by="rows", ascending=False)


def explore_table(conn, table_name):
    """Explore a specific table in detail."""
    try:
        # Get table structure
        columns = get_table_info(conn, table_name)

        # Get row count
        count = get_table_counts(conn, table_name)

        # Get sample data
        sample = get_sample_data(conn, table_name)

        # Get relationships
        relationships = find_related_tables(conn, table_name)

        # Print table information
        print(f"\n{'='*80}")
        print(f"TABLE: {table_name} ({count:,} rows)")
        print(f"{'='*80}")

        print("\nCOLUMNS:")
        print(tabulate(columns, headers="keys", tablefmt="psql"))

        if relationships:
            print("\nRELATIONSHIPS:")
            rel_df = pd.DataFrame(relationships)
            print(tabulate(rel_df, headers="keys", tablefmt="psql"))

        print("\nSAMPLE DATA:")
        print(tabulate(sample, headers="keys", tablefmt="psql", maxcolwidths=[None, 30, 30, 30, 30]))

    except Exception as e:
        print(f"Error exploring table {table_name}: {str(e)}")


def explorer_main(args):
    configs = check_and_download_chembl_db(version=args.version)
    with chembl_downloader.connect(version=configs["version"], prefix=configs["prefix"]) as conn:
        if args.list_tables:
            tables_df = list_all_tables(conn)
            print("\nTABLES IN CHEMBL DATABASE:")
            print(tabulate(tables_df, headers="keys", tablefmt="psql", showindex=False))

        elif args.table:
            explore_table(conn, args.table)

        elif args.search_column:
            results = search_tables_for_column(conn, args.search_column)
            if results:
                print(f"\nTables containing columns matching '{args.search_column}':")
                print(tabulate(pd.DataFrame(results), headers="keys", tablefmt="psql", showindex=False))
            else:
                print(f"No columns found matching '{args.search_column}'")

        elif args.query:
            try:
                result = pd.read_sql(args.query, conn)
                print("\nQuery Result:")
                print(tabulate(result, headers="keys", tablefmt="psql"))
            except Exception as e:
                print(f"Error executing query: {str(e)}")

        else:
            # Default: show database overview
            print("\nChEMBL DATABASE OVERVIEW")
            print("=" * 80)

            # Get tables and their counts
            tables_df = list_all_tables(conn)

            # Print top 20 tables by row count
            print("\nTop 20 tables by row count:")
            print(tabulate(tables_df.head(20), headers="keys", tablefmt="psql", showindex=False))

            # Suggest some common tables to explore
            print("\nSuggested tables to explore:")
            suggestions = [
                "molecule_dictionary - Main table for chemical compounds",
                "compound_structures - Chemical structures (SMILES, InChI)",
                "activities - Bioactivity measurements",
                "assays - Assay information",
                "target_dictionary - Biological targets",
                "docs - Document/publication information",
            ]

            for suggestion in suggestions:
                print(f"  • {suggestion}")

            print("\nTo explore a specific table:")
            print("  python chembl_schema_explorer.py --table <table_name>")
            print("\nTo list all tables:")
            print("  python chembl_schema_explorer.py --list-tables")
            print("\nTo search for columns:")
            print("  python chembl_schema_explorer.py --search-column <pattern>")

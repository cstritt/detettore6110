#!/usr/bin/env python3

import argparse
from find import main as find_main
from summarize import main as summarize_main

def main():
    parser = argparse.ArgumentParser(
        prog="detettore6110", 
        description="Infer insertion sequence polymorphisms and copy numbers from short-read sequencing data.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_find = subparsers.add_parser("find",
        help="Infer insertion sequence polymorphisms and copy numbers from short-read sequencing data.")
    p_find.set_defaults(func=find_main)

    p_sum = subparsers.add_parser("summarize", 
        help="Combine and summarize detettore6110 results.")
    p_sum.set_defaults(func=summarize_main)

    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return
    
    args.func()

if __name__ == "__main__":
    main()
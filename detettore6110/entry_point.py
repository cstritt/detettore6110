#!/usr/bin/env python3

import sys
from detettore6110.find import main as find_main
from detettore6110.summarize import main as summarize_main

def main():
    if len(sys.argv) < 2:
        print("usage: detettore6110 [-h] {find,summarize} ...")
        print("\nInfer insertion sequence polymorphisms and copy numbers from short-read sequencing data.")
        print("\npositional arguments:")
        print("  {find,summarize}")
        print("    find            Infer insertion sequence polymorphisms and copy numbers from short-read sequencing data.")
        print("    summarize       Combine and summarize detettore6110 results.")
        print("\noptions:")
        print("  -h, --help        show this help message and exit")
        return
    
    command = sys.argv[1]
    
    # Remove the command from sys.argv and let subcommand handle its own parsing
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    
    if command == "find":
        find_main()
    elif command == "summarize":
        summarize_main()
    else:
        print(f"Unknown command: {command}")
        print("Available commands: find, summarize")
        sys.exit(1)

if __name__ == "__main__":
    main()
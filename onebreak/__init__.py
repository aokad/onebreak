#! /usr/bin/env python

from .arg_parser import create_parser

def main():

    cparser = create_parser()
    args = cparser.parse_args()

    if vars(args) == {}:
        cparser.print_usage()
    else:
        args.func(args)


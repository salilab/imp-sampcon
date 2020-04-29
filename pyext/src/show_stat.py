from __future__ import print_function
from IMP import ArgumentParser

__doc__ = "Show fields in a PMI stat file."


def parse_args():
    parser = ArgumentParser(description="Show fields in a PMI stat file.")
    parser.add_argument("statfile", help="Name of the PMI stat file")
    return parser.parse_args()


def main():
    import IMP.pmi.output
    args = parse_args()

    p = IMP.pmi.output.ProcessOutput(args.statfile)
    fields = p.get_keys()
    print("\n".join(sorted(fields)))


if __name__ == '__main__':
    main()

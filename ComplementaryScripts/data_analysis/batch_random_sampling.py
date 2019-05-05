from random_sampling import run_random_sampling
import argparse

TIMEPOINT_DICT = {
    "M145": [21, 29, 33, 37, 41, 45, 49, 53, 57],
    "M1152":[33, 41, 45, 49, 53, 57, 61, 65]}


def run_batch(strains, iterations, processes):
    for s in strains:
        timepoints = TIMEPOINT_DICT[s]
        for t in timepoints:
            print(s, t, iterations, processes)
            run_random_sampling(s, t, iterations, processes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run random_sampling.py in batch")
    parser.add_argument("-s", "--strains", nargs = '+', required = True)
    parser.add_argument("-n", "--n_iter", help = "Number of iterations", type = int, required = True)
    parser.add_argument("-p", "--n_proc", help = "Number of processes", type = int, required = True)
    args = parser.parse_args()
    run_batch(args.strains, args.n_iter, args.n_proc)


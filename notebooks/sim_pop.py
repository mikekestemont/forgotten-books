import argparse
import itertools
import json
import multiprocessing as mp

from datetime import datetime

from tqdm import tqdm
import numpy as np
from numba import jit, njit, prange


class Parallel:
    def __init__(self, n_workers, n_tasks, pbar_position=0):
        self.pool = mp.Pool(n_workers)
        self._results = []
        self._pb = tqdm(total=n_tasks, position=pbar_position)

    def apply_async(self, fn, args=None):
        self.pool.apply_async(fn, args=args, callback=self._completed)

    def _completed(self, result):
        self._results.append(result)
        self._pb.update()

    def join(self):
        self.pool.close()
        self.pool.join()

    def result(self):
        self._pb.close()
        self.pool.close()
        return self._results

    
def reindex_array(x):
    _, x = np.unique(x, return_inverse=True)
    return x
    

@njit
def rand_choice_nb(arr, prob):
    return arr[np.searchsorted(np.cumsum(prob), np.random.random(), side="right")]


# @njit(parallel=True)
# def sample(population, adjacency_matrix):
#     new_population = np.zeros(population.shape[0])
#     for i in prange(adjacency_matrix.shape[0]):
#         neighborhood = population * adjacency_matrix[i]
#         counts = np.bincount(neighborhood[neighborhood > 0])
#         prob = counts / counts.sum()
#         traits = np.flatnonzero(counts)
#         new_population[i] = rand_choice_nb(traits, prob[prob > 0])
#     return new_population

@njit(parallel=True)
def sample(population, adjacency_matrix):
    new_population = np.zeros(population.shape[0])
    for i in prange(adjacency_matrix.shape[0]):
        neighborhood = population * adjacency_matrix[i]
        data = np.sort(neighborhood[neighborhood > 0])
        traits = np.unique(data)
        counts = np.zeros(traits.shape[0])
        count = 0
        j = 0
        for k in range(1, data.shape[0]):
            count += 1
            if data[k] > data[k - 1]:
                counts[j] = count
                count = 0
                j += 1
        counts[j] = count + 1
        prob = counts / counts.sum()
        new_population[i] = rand_choice_nb(traits, prob)
    return new_population


def generate_network(n, k):
    """Returns a graph with n nodes and k nearest neighbors in a ring topology."""
    adjacency_matrix = np.ones((n, n), dtype=np.int8)
    if n == k:
        np.fill_diagonal(adjacency_matrix, 0)
        return adjacency_matrix

    adjacency_matrix = adjacency_matrix * 0
    nodes = np.arange(n)
    for i in range(1, k // 2 + 1):
        targets = np.take(nodes, range(i, i + n), mode="wrap")
        adjacency_matrix[nodes, targets] = 1
        adjacency_matrix[targets, nodes] = 1
    return adjacency_matrix


class NeutralModel:
    """
    Agent-based Evolutionary Model.
    """

    def __init__(
        self, n=100_000, mu=0.0005, init_num_traits=2, seed=None, disable_pbar=False
    ):
        self.n = n
        self.mu = mu
        self.init_num_traits = init_num_traits
        self.disable_pbar = disable_pbar

        # random number generator
        self.rnd = np.random.RandomState(seed)

        # initialize population:
        self.population = np.ones(self.n, dtype=np.int64)
        for i in range(1, self.init_num_traits + 1):
            self.population[
                int((i - 1) * self.n / self.init_num_traits) : int(
                    i * self.n / self.init_num_traits
                )
            ] = i

        self.num_traits = len(np.unique(self.population)) + 1
        assert self.num_traits == self.init_num_traits + 1

    def __repr__(self):
        return f"agents={self.n}, mu={self.mu}"

    def fit(self):
        """
        Burnin
        """
        iterations = 0
        with tqdm(desc=f"burn-in: {self}", disable=self.disable_pbar, position=1) as pb:
            while self.population.min() <= self.init_num_traits:
                iterations += 1
                self.step()
                pb.update()
        # print(f'burn-in halted after {iterations} iterations')
        self.population = reindex_array(self.population)
        return self

    def _sample(self):
        traits, counts = np.unique(self.population, return_counts=True)
        return self.rnd.choice(traits, self.n, replace=True, p=counts / counts.sum())

    def step(self):
        """
        One forward pass
        """
        self.population = self._sample()

        # Check who was in reality innovating
        innovating_population = self.rnd.rand(self.n) < self.mu
        n_innovations = innovating_population.sum()

        # Assign new traits to the innovating individuals
        self.population[innovating_population] = np.arange(
            self.num_traits, self.num_traits + n_innovations
        )

        # Update the total number of traits
        self.num_traits += n_innovations

    def postfit(self, iterations=10_000):
        for _ in tqdm(range(iterations), desc=f"post-burn-in: {self}", disable=self.disable_pbar, position=2):
            self.step()
        return self


class NetworkModel(NeutralModel):
    """
    Agent-based Evolutionary Model.
    """

    def __init__(
        self,
        n=100_000,
        mu=0.0005,
        init_num_traits=2,
        k=3,
        disable_pbar=False,
        seed=None,
    ):
        super().__init__(
            n=n,
            mu=mu,
            init_num_traits=init_num_traits,
            disable_pbar=disable_pbar,
            seed=seed,
        )
        self.k = k

        # initialize network topology
        self.adjacency_matrix = generate_network(self.n, self.k)
        assert self.adjacency_matrix.shape == (self.n, self.n)

    def __repr__(self):
        return f"agents={self.n}, mu={self.mu}, k={self.k}"

    def _sample(self):
        return sample(self.population, self.adjacency_matrix).astype(np.int64)


def simulate(agents, k, mu, postfit_iterations):
    if agents == k:
        simulator = NeutralModel(n=agents, mu=mu, disable_pbar=True)
    else:
        simulator = NetworkModel(n=agents, mu=mu, k=k, disable_pbar=False)
    simulator.fit()
    simulator.postfit(postfit_iterations)
    return np.array([agents, k, mu]), simulator.population


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--agents", type=int, nargs="+", default=(1_000,))
    parser.add_argument("-k", "--degree", type=int, nargs="+", default=(2, 10, 100))
    parser.add_argument("-s", "--simulations", type=int, default=1_000)
    parser.add_argument("-r", "--reference", action="store_true")
    parser.add_argument("-m", "--mu", type=float, nargs="+", default=None)
    parser.add_argument("-i", "--iterations", type=int, default=10_000)
    parser.add_argument("-w", "--workers", type=int, default=1)
    args = parser.parse_args()

    if args.mu is None:
        mu = [0.005] * len(args.agents) #tuple([50 / agents for agents in args.agents])
        args.__dict__["mu"] = mu

    pool = Parallel(
        args.workers, args.simulations * len(args.agents) * len(args.degree),
        pbar_position=0
    )
    for i in range(args.simulations):
        if not args.reference:
            for j, agents in enumerate(args.agents):
                mu = args.mu[j]
                for degree in args.degree:
                    pool.apply_async(simulate, args=(agents, degree, mu, args.iterations))
        else:
            for j, agents in enumerate(args.agents):
                mu = args.mu[j]
                pool.apply_async(simulate, args=(agents, agents, mu, args.iterations))
    pool.join()

    params, populations = zip(*pool.result())
    max_len = max(len(population) for population in populations)
    populations = [
        np.pad(population, (0, max_len - len(population)), "constant")
        for population in populations
    ]
    params, populations = np.vstack(params), np.vstack(populations)

    now = datetime.now().strftime("%Y%m%d%H%M%S")
    np.savez_compressed(
        f"{now}.npz",
        params=params.astype(np.float64),
        populations=populations,
    )

    with open(f"{now}.params.json", "w") as f:
        json.dump(args.__dict__, f)

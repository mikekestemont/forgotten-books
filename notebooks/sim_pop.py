import argparse
import json
import itertools

from tqdm import tqdm
import numpy as np
from numba import jit, njit, prange


@njit
def rand_choice_nb(arr, prob):
    return arr[np.searchsorted(np.cumsum(prob), np.random.random(), side="right")]


@njit(parallel=True)
def sample(population, adjacency_matrix):
    new_population = np.zeros(population.shape[0])
    for i in prange(adjacency_matrix.shape[0]):
        neighborhood = population * adjacency_matrix[i]
        counts = np.bincount(neighborhood[neighborhood > 0])
        prob = counts / counts.sum()
        traits = np.flatnonzero(counts)
        new_population[i] = rand_choice_nb(traits, prob[prob > 0])
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

    def __init__(self, n=100_000, mu=0.0005, init_num_traits=2, seed=None):
        self.n = n
        self.mu = mu
        self.init_num_traits = init_num_traits

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

    def fit(self):
        """
        Burnin
        """
        iterations = 0
        with tqdm(desc="burn-in") as pb:
            while self.population.min() <= self.init_num_traits:
                iterations += 1
                self.step()
                pb.update()
        # print(f'burn-in halted after {iterations} iterations')
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
        for _ in tqdm(range(iterations), desc="post-burn-in"):
            self.step()
        return self


class NetworkModel(NeutralModel):
    """
    Agent-based Evolutionary Model.
    """

    def __init__(self, n=100_000, mu=0.0005, init_num_traits=2, k=3, seed=None):
        super().__init__(n=n, mu=mu, init_num_traits=init_num_traits, seed=seed)
        self.k = k

        if self.n == self.k:
            return NeutralModel(n, mu, init_num_traits, seed=seed)

        # initialize network topology
        self.adjacency_matrix = generate_network(self.n, self.k)
        assert self.adjacency_matrix.shape == (self.n, self.n)

    def _sample(self):
        return sample(self.population, self.adjacency_matrix).astype(np.int64)


# parameter ranges:
num_agents = (10000, )
ks = (10, 25)
num_simulations = 5

# fixed parameters:
num_postfit = 100
init_num_traits = 2
mu = 0.0005
outfile = 'populations.json'

# first, network simulations:
populations = {k:[] for k in ks}
network_combs = list(itertools.product(*(num_agents, ks,
                                list(range(num_simulations)))))
for idx, (agents, k, sim) in enumerate(network_combs):
    print(f'# agents={agents} | k={k} | simulation {sim + 1}/{num_simulations} | network-exp {idx + 1}/{len(network_combs)}')
    p = NetworkModel(n=agents, k=k, mu=mu, init_num_traits=init_num_traits)\
                    .fit().postfit(num_postfit).population
    populations[k].append(p.tolist())

# secondly, neutral simulations:
reference_combs = list(itertools.product(*(num_agents,
                                list(range(num_simulations)))))
populations['neutral'] = []
for idx, (na, sim) in enumerate(reference_combs):
    print(f'# agents={agents} | k=neutral | simulation {sim + 1}/{num_simulations} | neutral-exp {idx + 1}/{len(reference_combs)}')
    p = NeutralModel(n=agents, mu=mu, init_num_traits=init_num_traits)\
                    .fit().postfit(num_postfit).population
    populations['neutral'].append(p.tolist())

with open(outfile, 'w') as f:
    f.write(json.dumps(populations))
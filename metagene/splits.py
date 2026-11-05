"""Observation splits that preserve the scientific unit of generalization."""

from dataclasses import dataclass
import numpy as np


@dataclass
class Split:
    train: np.ndarray
    val: np.ndarray
    test: np.ndarray
    name: str
    seed: int
    discarded: np.ndarray | None = None


def _grouped_partition(pairs, indices, val_fraction, test_fraction, rng):
    canonical = np.sort(np.asarray(pairs)[indices], axis=1)
    _, group = np.unique(canonical, axis=0, return_inverse=True)
    groups = rng.permutation(group.max() + 1)
    n_test = round(len(groups) * test_fraction)
    n_val = round(len(groups) * val_fraction)
    test = indices[np.isin(group, groups[:n_test])]
    val = indices[np.isin(group, groups[n_test:n_test + n_val])]
    train = indices[np.isin(group, groups[n_test + n_val:])]
    return train, val, test


def pair_random_split(pairs, seed=0, val_fraction=0.15, test_fraction=0.15):
    """Keep all observations of an unordered gene pair in the same partition."""
    indices = np.arange(len(pairs))
    train, val, test = _grouped_partition(
        pairs, indices, val_fraction, test_fraction, np.random.default_rng(seed)
    )
    return Split(train, val, test, "pair_random", seed)


def gene_disjoint_split(pairs, seed=0, val_fraction=0.2, test_fraction=0.2):
    """Both test genes are absent from training; mixed-partition pairs are unused.

    Fractions allocate genes, not observations. The retained pair counts are
    consequently smaller than a random pair split and must be reported.
    """
    pairs = np.asarray(pairs)
    genes = np.random.default_rng(seed).permutation(np.unique(pairs))
    n_test = round(len(genes) * test_fraction)
    n_val = round(len(genes) * val_fraction)
    test_genes = genes[:n_test]
    val_genes = genes[n_test:n_test + n_val]
    train_genes = genes[n_test + n_val:]
    train_mask = np.isin(pairs, train_genes).all(axis=1)
    val_mask = np.isin(pairs, val_genes).all(axis=1)
    test_mask = np.isin(pairs, test_genes).all(axis=1)
    return Split(
        np.flatnonzero(train_mask), np.flatnonzero(val_mask),
        np.flatnonzero(test_mask), "gene_disjoint", seed,
        np.flatnonzero(~(train_mask | val_mask | test_mask)),
    )


def leave_context_out_split(pairs, context_ids, held_out_context, seed=0,
                            val_fraction=0.15, validation_context=None):
    """Hold out a complete context; split remaining pairs for model selection.

    A named validation context instead permits selection across contexts.
    The same pair may occur in training and held-out contexts by design.
    """
    contexts = np.asarray(context_ids)
    test = np.flatnonzero(contexts == held_out_context)
    remaining = np.flatnonzero(contexts != held_out_context)
    if validation_context is None:
        train, val, _ = _grouped_partition(
            pairs, remaining, val_fraction, 0, np.random.default_rng(seed)
        )
    else:
        val = remaining[contexts[remaining] == validation_context]
        train = remaining[contexts[remaining] != validation_context]
    return Split(train, val, test, f"leave_context_{held_out_context}", seed)


def restrict_split(split, indices):
    """Restrict every partition to an explicitly selected measured dataset."""
    indices = np.asarray(indices)
    return Split(*(x[np.isin(x, indices)] for x in (split.train, split.val, split.test)),
                 split.name, split.seed, split.discarded)

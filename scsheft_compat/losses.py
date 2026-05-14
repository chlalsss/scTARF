"""Compatibility losses for CSUBioGroup/scSHEFT.

The public repository imports ``losses.py`` from ``scSHEFT.py``, but the cloned
repository may not include that file. These implementations provide the symbols
used by the official training loop without modifying the upstream checkout.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class InfoNCE(nn.Module):
    def __init__(self, batch_size: int, tau: float = 0.4):
        super().__init__()
        self.batch_size = batch_size
        self.tau = tau

    def forward(self, features: torch.Tensor, positive_features: torch.Tensor) -> torch.Tensor:
        z1 = F.normalize(features, p=2, dim=1)
        z2 = F.normalize(positive_features, p=2, dim=1)
        logits = torch.matmul(z1, z2.t()) / self.tau
        labels = torch.arange(logits.shape[0], device=logits.device)
        return F.cross_entropy(logits, labels)


class L1Regularization(nn.Module):
    def __init__(self, weight: float = 0.0):
        super().__init__()
        self.weight = weight

    def forward(self, model: nn.Module) -> torch.Tensor:
        if self.weight == 0:
            try:
                param = next(model.parameters())
                return param.new_tensor(0.0)
            except StopIteration:
                return torch.tensor(0.0)

        loss = None
        for param in model.parameters():
            term = param.abs().sum()
            loss = term if loss is None else loss + term
        if loss is None:
            return torch.tensor(0.0)
        return self.weight * loss


def zero_center(features: torch.Tensor) -> torch.Tensor:
    return features.mean(dim=0).pow(2).mean()


def non_corr(features: torch.Tensor) -> torch.Tensor:
    if features.shape[0] <= 1:
        return features.new_tensor(0.0)
    z = features - features.mean(dim=0, keepdim=True)
    z = z / (z.std(dim=0, keepdim=True) + 1e-6)
    corr = torch.matmul(z.t(), z) / max(features.shape[0] - 1, 1)
    eye = torch.eye(corr.shape[0], device=corr.device, dtype=corr.dtype)
    return ((corr - eye) ** 2).mean()


def max_var(features: torch.Tensor) -> torch.Tensor:
    std = torch.sqrt(features.var(dim=0) + 1e-6)
    return F.relu(1.0 - std).mean()

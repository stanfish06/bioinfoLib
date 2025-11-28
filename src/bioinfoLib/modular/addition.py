# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import time

import numpy as np
from max import engine
from max.driver import CPU
from max.dtype import DType
from max.graph import DeviceRef, Graph, TensorType, ops


def add_tensors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    input_type = TensorType(DType.float32, shape=a.shape, device=DeviceRef.CPU())
    with Graph("add", input_types=(input_type, input_type)) as graph:
        lhs, rhs = graph.inputs
        out = ops.add(lhs, rhs)
        graph.output(out)
    session = engine.InferenceSession(devices=[CPU()])
    model = session.load(graph)
    output = model.execute(a, b)[0]
    return output.to_numpy()


def main():
    a = np.random.randn(100).reshape([20, -1]).astype(np.float32)
    b = np.random.randn(100).reshape([20, -1]).astype(np.float32)
    start = time.time()
    add_tensors(a, b)
    end = time.time()
    elapsed = end - start
    print(f"Elapsed: {elapsed:.4f} seconds")

    start = time.time()
    a + b
    end = time.time()
    elapsed = end - start
    print(f"Elapsed: {elapsed:.4f} seconds")


if __name__ == "__main__":
    main()

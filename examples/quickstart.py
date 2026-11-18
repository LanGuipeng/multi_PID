from itertools import product

from mar_certified import Distribution, query

counts = [8, 12, 32, 48, 76, 4, 19, 1]
distribution = Distribution(dict(zip(product((0, 1), repeat=3), counts)), kind="counts")
result = query(distribution, [[1], [2]], tol="1e-12")
assert result["tolerance_met"]
print(f"Redundancy: {result['approximate']:.15f} bit")
print(f"Certified width: {float(result['interval'].width):.3e} bit")

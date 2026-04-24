from collections import defaultdict


def make_vector(R):
    """Factory for a module over semiring R.

    The returned `Vector` class is NOT itself a semiring — it's a vector
    space / module: Vector + Vector, Vector * scalar(R), Vector @ Vector,
    Vector / scalar(R). Elements are `R.chart()` maps.

    Requires R to provide `chart()`, `lift(x)`, `zero`. For `to_real()`,
    values in the vector must support `.to_real()`. For `is_zero()`, values
    must support `.is_zero()`.
    """

    class Vector:
        def __init__(self, vals=None):
            self.x = R.chart()
            if vals is not None:
                self.x.update(vals)

        def __mul__(self, x):
            assert isinstance(x, R), x
            c = self.__class__()
            for k, v in self.x.items():
                c.x[k] = v * x
            return c

        __rmul__ = __mul__

        def __truediv__(self, x):
            assert isinstance(x, R), x
            c = self.__class__()
            for k, v in self.x.items():
                c.x[k] = v / x
            return c

        __div__ = __truediv__

        def __sub__(self, x):
            return self + R.lift(-1) * x

        def __matmul__(self, w):
            v = R.zero
            for k in self:
                v += self[k] * w[k]
            return v

        def __add__(self, x):
            assert isinstance(x, Vector), x
            c = self.__class__()
            for k, v in self.x.items():
                c.x[k] += v
            for k, v in x.x.items():
                c.x[k] += v
            return c

        def is_zero(self):
            return all(y.is_zero() for y in self.x.values())

        def __repr__(self):
            return f'Vector({type(R).__name__}, {dict(self.x)})'

        def __str__(self):
            innards = '\n'.join('  %r: %r,' % (k, v) for k, v in sorted(self.x.items()))
            if innards:
                innards = '{\n%s\n}' % innards
            return f'Vector({type(R).__name__}, {innards})'

        def __getitem__(self, key):
            return self.x[key]

        def to_real(self):
            x = defaultdict(float)
            for k, v in self.x.items():
                x[k] = v.to_real()
            return x

        def __iter__(self):
            return iter(self.x)

    Vector.zero = Vector()

    return Vector

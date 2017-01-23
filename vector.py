from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 30


class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel compoenent'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = 'No unique orhtogonal compoenent'
    ONLY_DEFINED_IN_TWO_THREE_DIM_MSG = 'Only defined in two three dim vectors'

    def __init__(self, coords):

        try:
            if not coords:
                raise ValueError
            self.coords = tuple([Decimal(x) for x in coords])
            self.dimension = len(coords)

        except ValueError:
            raise ValueError('The coords must be nonempty')

        except TypeError:
            raise TypeError('The coords must be an iterable')

    def __getitem__(self, i):
        '''get the coordinate at specified position'''
        return self.coords[i]

    def __str__(self):
        '''print out vector in a nice way'''
        return 'Vector: {}'.format(self.coords)

    def __eq__(self, v):
        '''two vectors are equal if their coordinates are equal'''
        return self.coords == v.coords

    def plus(self, v):
        '''add a Vector to a Vector'''
        new_coords = [x + y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def minus(self, v):
        '''subtract a vector from a vector'''
        new_coords = [x - y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def times_scalar(self, scalar):
        '''Scalar Multiplication'''
        new_coords = [Decimal(scalar) * x for x in self.coords]
        return Vector(new_coords)

    def magnitude(self):
        '''calculate length of a vector (pythagoras)'''
        coords_squared = [x**2 for x in self.coords]
        return Decimal(sqrt(sum(coords_squared)))

    def normalized(self):
        '''calculate vector of magnitude 1, pointing in the same direction'''
        try:
            magnitude = Decimal(self.magnitude())
            return self.times_scalar(Decimal('1.0') / magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        '''calculate Dot-Product (Inner-/Scalar-Product) of two vectors'''
        return sum([x * y for x, y in zip(self.coords, v.coords)])

    def angle_with(self, v, in_degrees=False):
        '''calculate angle between two vectors'''
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            dot_product = round(u1.dot(u2), 6)
            angle_in_radians = acos(dot_product)

            if in_degrees:
                degrees_per_radian = 180. / pi
                return Decimal(angle_in_radians * degrees_per_radian)
            else:
                return Decimal(angle_in_radians)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')

    def is_parallel_to(self, v):
        '''either one is zero vector or they are pointing in the same/ opposite
        direction'''
        return (self.is_zero() or
                v.is_zero() or
                self.angle_with(v) == 0 or
                self.angle_with(v) == pi)

    def is_orthogonal_to(self, v, tolerance=1e-10):
        '''two vectors are orthogonal if their Dot-Product is 0'''
        return abs(self.dot(v)) < tolerance

    def is_zero(self, tolerance=1e-10):
        '''test if magnitude is smaller than tolerance'''
        return self.magnitude() < tolerance

    def component_parallel_to(self, basis):
        '''calculate the projection of the vector onto a basis'''
        try:
            unit_vector = basis.normalized()
            weight = self.dot(unit_vector)
            return unit_vector.times_scalar(weight)
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def component_orthogonal_to(self, basis):
        '''calculate component that is orthogonal to the basis'''
        try:
            basis_projection = self.component_parallel_to(basis)
            return self.minus(basis_projection)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def cross_product_with(self, v):
        '''find vector which is orthogonal to both vectors'''
        try:
            x1, y1, z1 = self.coords
            x2, y2, z2 = v.coords
            new_coords = [y1 * z2 - y2 * z1,
                          -(x1 * z2 - x2 * z1),
                          x1 * y2 - x2 * y1]
            return Vector(new_coords)
        except ValueError as e:
            msg = str(e)
            # if msg == 'need more than 2 values to unpack':
            if msg == 'not enough values to unpack (expected 3, got 2)':
                self_embedded_in_R3 = Vector(self.coords + ('0',))
                v_embedded_in_R3 = Vector(v.coords + ('0',))
                return self_embedded_in_R3.cross_product_with(v_embedded_in_R3)
            elif (msg == 'too many values to unpack' or
                    msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
            else:
                raise e

    def area_of_parallelogram_with(self, v):
        '''Area of parallelogram == magnitude of the Cross-Product'''
        cross_product = self.cross_product_with(v)
        return cross_product.magnitude()

    def area_of_triangle_with(self, v):
        '''area of triangle == half magnitude of the Cross-Product'''
        return Decimal('0.5') * self.area_of_parallelogram_with(v)


if __name__ == '__main__':

    print('##################################')
    print('Quiz: Plus, Minus, Scalar Multiply')

    v1 = Vector([8.218, -9.341])
    v2 = Vector([-1.129, 2.111])
    result = v1.plus(v2)
    print(result)

    v1 = Vector([7.119, 8.215])
    v2 = Vector([-8.223, 0.878])
    result = v1.minus(v2)
    print(result)

    v1 = Vector([1.671, -1.012, -0.318])
    s1 = 7.41
    result = v1.times_scalar(s1)
    print(result)

    print('##################################')
    print('Quiz: Coding Magnitude & Direction')

    v1 = Vector([-0.221, 7.437])
    result = v1.magnitude()
    print(result)

    v1 = Vector([8.813, -1.331, -6.247])
    result = v1.magnitude()
    print(result)

    v1 = Vector([5.581, -2.136])
    result = v1.normalized()
    print(result)

    v1 = Vector([1.996, 3.108, -4.554])
    result = v1.normalized()
    print(result)

    print('################################')
    print('Quiz: Coding Dot Product & Angle')

    v1 = Vector([7.887, 4.138])
    w1 = Vector([-8.802, 6.776])
    result = v1.dot(w1)
    print(result)

    v2 = Vector([-5.955, -4.904, -1.874])
    w2 = Vector([-4.496, -8.755, 7.103])
    result = v2.dot(w2)
    print(result)

    v3 = Vector([3.183, -7.627])
    w3 = Vector([-2.668, 5.319])
    result = v3.angle_with(w3)
    print(result)

    v4 = Vector([7.35, 0.221, 5.188])
    w4 = Vector([2.751, 8.259, 3.985])
    result = v4.angle_with(w4, True)
    print(result)

    print('###################################')
    print('Quiz: Checking Parallel, Orthogonal')

    v1 = Vector([-7.579, -7.88])
    w1 = Vector([22.737, 23.64])
    print('parallel: {0}'.format(v1.is_parallel_to(w1)))
    print('orthogonal: {0}'.format(v1.is_orthogonal_to(w1)))

    v2 = Vector([-2.029, 9.97, 4.172])
    w2 = Vector([-9.231, -6.639, -7.245])
    print('parallel: {0}'.format(v2.is_parallel_to(w2)))
    print('orthogonal: {0}'.format(v2.is_orthogonal_to(w2)))

    v3 = Vector([-2.328, -7.284, -1.214])
    w3 = Vector([-1.821, 1.072, -2.94])
    print('parallel: {0}'.format(v3.is_parallel_to(w3)))
    print('orthogonal: {0}'.format(v3.is_orthogonal_to(w3)))

    v4 = Vector([2.118, 4.827])
    w4 = Vector([0, 0])
    print('parallel: {0}'.format(v4.is_parallel_to(w4)))
    print('orthogonal: {0}'.format(v4.is_orthogonal_to(w4)))

    print('###############################')
    print('Quiz: Coding Vector Projections')

    v1 = Vector([3.039, 1.879])
    b1 = Vector([0.825, 2.036])
    print(v1.component_parallel_to(b1))

    v2 = Vector([-9.88, -3.264, -8.159])
    b2 = Vector([-2.155, -9.353, -9.473])
    print(v2.component_orthogonal_to(b2))

    v3 = Vector([3.009, -6.172, 3.692, -2.51])
    b3 = Vector([6.404, -9.144, 2.759, 8.718])
    print(v3.component_parallel_to(b3))
    print(v3.component_orthogonal_to(b3))

    print('###########################')
    print('Quiz: Coding Cross Products')

    v1 = Vector([8.462, 7.893, -8.187])
    w1 = Vector([6.984, -5.975, 4.778])
    print(v1.cross_product_with(w1))

    v2 = Vector([-8.987, -9.838, 5.031])
    w2 = Vector([-4.268, -1.861, -8.866])
    print(v2.area_of_parallelogram_with(w2))

    v3 = Vector([1.5, 9.547, 3.691])
    w3 = Vector([-6.007, 0.124, 5.772])
    print(v3.area_of_triangle_with(w3))

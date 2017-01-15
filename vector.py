from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 3


class Vector(object):
    '''Defines a Vector with its properties and methods'''

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

    def __str__(self):
        return 'Vector: {}'.format(self.coords)

    def __eq__(self, v):
        return self.coords == v.coords

    def plus(self, v):
        '''add a Vector to a Vector

        Args:
            Vector Object which should be added
        Returns:
            Vector Object
        '''
        new_coords = [x+y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def minus(self, v):
        '''subtract a vector from a vector

        Args:
            Vector Object which should be subtracted
        Returns:
            Vector Object
        '''
        new_coords = [x-y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def times_scalar(self, scalar):
        '''Scalar Multiplication

        Args:
            scalar value
        Returns:
            Vector Object
        '''
        new_coords = [Decimal(scalar)*x for x in self.coords]
        return Vector(new_coords)

    def magnitude(self):
        '''calculate length of a vector (pythagoras)

        magnitude of a vector is the distance between start and end
        coords

        Returns:
            scalar value
        '''
        coords_squared = [x**2 for x in self.coords]
        return Decimal(sqrt(sum(coords_squared)))

    def normalized(self):
        '''finding the unit vector pointing in the same direction

        A unit Vector is a vector whose magnitude is 1.

        Args:
            Vector Object itself
        Returns:
            Vector Object
        '''
        try:
            magnitude = Decimal(self.magnitude())
            return self.times_scalar(Decimal('1.0')/magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        '''calculate the dot product of two vectors

        Alternative names are Inner Product or Scalar Product

        Args:
            v: Vector Object
        Returns:
            scalar value <class 'decimal.Decimal'>
        '''
        return sum([x*y for x, y in zip(self.coords, v.coords)])

    def angle_with(self, v, in_degrees=False):
        '''calculate angle between two vectors

        Args:
            v: Vector Object
            in_degrees: boolean value. false(default)=return angle in randians,
                true=return angle in degrees
        Returns:
            scalar value <class 'decimal.Decimal'>
        '''
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            dot_product = round(u1.dot(u2), 3)
            angle_in_radians = acos(dot_product)

            if in_degrees:
                degrees_per_radian = 180./pi
                return Decimal(angle_in_radians*degrees_per_radian)
            else:
                return Decimal(angle_in_radians)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')

    def is_parallel_to(self, v):
        '''check if two vectors are parallel

        two vectors are parallel if one vector is the zero vector. However if
        they are both non zero, they both must be pointing in the same
        direction or in completly opposite direction

        args:
            v: Vector Object
        returns:
            boolean value
        '''
        return (self.is_zero() or
                v.is_zero() or
                self.angle_with(v) == 0 or
                self.angle_with(v) == pi)

    def is_orthogonal_to(self, v, tolerance=1e-10):
        '''check if two vectors are orthogonal

        two vectors are orthogonal if their dot Product is 0

        args:
            v: Vector Object
            tolerance: float value to balance rounding issues
        returns:
            boolean value <class 'bool'>
        '''
        return abs(self.dot(v)) < tolerance

    def is_zero(self, tolerance=1e-10):
        '''test if magnitude is smaller than tolerance

        args:
            tolerance: float value, default is 1e-10
        returns:
            boolean value <class 'bool'>
        '''
        return self.magnitude() < tolerance

    def component_parallel_to(self, basis):
        '''calculate the projection of the vector onto a basis

        args:
            basis: Vector Object
        returns:
            Vector Object
        '''
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
        '''calculate component of the vector that is orthogonal to the basis

        calculate the difference between the the vector and its projection
        onto the basis

        args:
            basis: Vector Object
        returns:
            Vector Object
        '''
        try:
            basis_projection = self.component_parallel_to(basis)
            return self.minus(basis_projection)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def cross_product_with(self, v):
        '''calculate the cross product of two vectors

        args:
            v: Vector Object
        returns:
            Vector Object
        '''
        try:
            x1, y1, z1 = self.coords
            x2, y2, z2 = v.coords
            new_coords = [y1*z2-y2*z1, -(x1*z2-x2*z1), x1*y2-x2*y1]
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
        '''calculate the area of parallelogram spanned by two vectors

        args:
            v: Vector Object
        returns:
            scalar value
        '''
        cross_product = self.cross_product_with(v)
        return cross_product.magnitude()

    def area_of_triangle_with(self, v):
        ''' calculate the area between two vectors

        half size of the area of parallelogram spanned by two vectors

        args:
            v: Vector Object
        returns:
            scalar Value
        '''
        return Decimal('0.5') * self.area_of_parallelogram_with(v)

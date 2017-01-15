from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class Plane(object):
    '''Vector representation of a Plane in three dimension'''

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        '''Create a Plane Object

        consists of a normal vector, constant term and basepoint
        Ax + By + Cz = k
        Vector([A, B, C]) represents a normal Vector

        Args:
            normal_vector: Vector Object
            constant_term:
        '''
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def set_basepoint(self):
        '''calculate basepoint

        find the the first non zero coordinate, either (0, k/B) or (k/A, 0)
        '''
        try:
            n = self.normal_vector.coords
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector.coords

        try:
            initial_index = Plane.first_nonzero_index(n)
            terms = []
            for i in range(self.dimension):
                if round(n[i], num_decimal_places) != 0:
                    terms.append(write_coefficient(
                        n[i],
                        is_initial_term=(i == initial_index))
                        + 'x_{}'.format(i+1))
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel_to(self, plane2):
        '''determine if this plane is parallel to plane2

        two planes are parallel if their normal vectors are parallel

        args:
            plane2: Plane Object
        '''
        return self.normal_vector.is_parallel_to(plane2.normal_vector)

    def __eq__(self, plane2):
        '''determine if this plane is equal to plane2

        two planes are equal, if the vector connecting one point on each plane
        is orthogonal to the planes normal vectors

        args:
            plane2: Plane Object
        '''
        # normal vector have to be parallel in order to be equal
        if not self.is_parallel_to(plane2):
            return False

        v_connect = self.basepoint.minus(plane2.basepoint)
        return v_connect.is_orthogonal_to(self.normal_vector)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

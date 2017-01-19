from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class Plane(object):
    '''Vector representation of a Plane in three dimension'''

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        '''Create a Plane Object

        Ax + By + Cz = k
        Vector([A, B, C]) represents a normal Vector'''
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def __getitem__(self, i):
        return self.normal_vector[i]

    def set_basepoint(self):
        '''find the first non zero coordinate, either (0, k/B) or (k/A, 0)'''
        try:
            n = self.normal_vector
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

        n = self.normal_vector

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
        '''two planes are parallel if their normal vectors are parallel'''
        return self.normal_vector.is_parallel_to(plane2.normal_vector)

    def __eq__(self, plane2):
        '''two planes are equal, if the vector connecting one point on each
        plane is orthogonal to the planes normal vectors'''
        # normal vector have to be parallel in order to be equal
        # if not self.is_parallel_to(plane2):
        #     return False

        # v_connect = self.basepoint.minus(plane2.basepoint)
        # return v_connect.is_orthogonal_to(self.normal_vector)
        if self.normal_vector.is_zero():
            if not plane2.normal_vector.is_zero():
                return False

            diff = self.constant_term - plane2.constant_term
            return MyDecimal(diff).is_near_zero()

        elif plane2.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(plane2):
            return False

        basepoint_difference = self.basepoint.minus(plane2.basepoint)
        return basepoint_difference.is_orthogonal_to(self.normal_vector)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

if __name__ == '__main__':

    print('################################')
    print('Quiz: Planes in 3 Dimensions - 2')

    p1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
    p2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
    print('equal: {0}'.format(p1 == p2))
    print('parallel but unequal: {0}'.format(p1.is_parallel_to(p2)))

    p3 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
    p4 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)
    print('equal: {0}'.format(p3 == p4))
    print('parallel but unequal: {0}'.format(p3.is_parallel_to(p4)))

    p5 = Plane(Vector([-7.926, 8.625, -7.212]), -7.952)
    p6 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)
    print('equal: {0}'.format(p5 == p6))
    print('parallel but unequal: {0}'.format(p5.is_parallel_to(p6)))

from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class Hyperplane(object):
    '''Vector representation of a Hyperplane'''

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    EITHER_DIM_OR_NORMAL_VEC_MUST_BE_PROVIDED_MSG = 'Either the dimension or the normal vector must be provided'

    def __init__(self, dimension=None, normal_vector=None, constant_term=None):
        '''Create a Hyperplane Object'''
        if not dimension and not normal_vector:
            raise Exception(self.EITHER_DIM_OR_NORMAL_VEC_MUST_BE_PROVIDED_MSG)

        elif not normal_vector:
            self.dimension = dimension
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        else:
            self.dimension = normal_vector.dimension
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def __getitem__(self, i):
        return self.normal_vector[i]

    def set_basepoint(self):
        '''find the first non zero coordinate'''
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension

            initial_index = Hyperplane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c / initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Hyperplane.NO_NONZERO_ELTS_FOUND_MSG:
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
            initial_index = Hyperplane.first_nonzero_index(n)
            terms = []
            for i in range(self.dimension):
                if round(n[i], num_decimal_places) != 0:
                    terms.append(write_coefficient(
                        n[i],
                        is_initial_term=(i == initial_index)) +
                        'x_{}'.format(i + 1))
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
        raise Exception(Hyperplane.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel_to(self, Hyperplane2):
        '''two Hyperplanes are parallel if their normal vectors are parallel'''
        return self.normal_vector.is_parallel_to(Hyperplane2.normal_vector)

    def __eq__(self, Hyperplane2):
        '''two Hyperplanes are equal, if the vector connecting one point on each
        Hyperplane is orthogonal to the Hyperplanes normal vectors'''
        # normal vector have to be parallel in order to be equal
        # if not self.is_parallel_to(Hyperplane2):
        #     return False

        # v_connect = self.basepoint.minus(Hyperplane2.basepoint)
        # return v_connect.is_orthogonal_to(self.normal_vector)
        if self.normal_vector.is_zero():
            if not Hyperplane2.normal_vector.is_zero():
                return False

            diff = self.constant_term - Hyperplane2.constant_term
            return MyDecimal(diff).is_near_zero()

        elif Hyperplane2.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(Hyperplane2):
            return False

        basepoint_difference = self.basepoint.minus(Hyperplane2.basepoint)
        return basepoint_difference.is_orthogonal_to(self.normal_vector)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


if __name__ == '__main__':

    print('################################')
    print('Quiz: Hyperplanes in 3 Dimensions - 2')

    p1 = Hyperplane(normal_vector=Vector([-0.412, 3.806, 0.728]), constant_term=-3.46)
    p2 = Hyperplane(normal_vector=Vector([1.03, -9.515, -1.82]), constant_term=8.65)
    print('equal: {0}'.format(p1 == p2))
    print('parallel but unequal: {0}'.format(p1.is_parallel_to(p2)))

    p3 = Hyperplane(normal_vector=Vector([2.611, 5.528, 0.283]), constant_term=4.6)
    p4 = Hyperplane(normal_vector=Vector([7.715, 8.306, 5.342]), constant_term=3.76)
    print('equal: {0}'.format(p3 == p4))
    print('parallel but unequal: {0}'.format(p3.is_parallel_to(p4)))

    p5 = Hyperplane(normal_vector=Vector([-7.926, 8.625, -7.212]), constant_term=-7.952)
    p6 = Hyperplane(normal_vector=Vector([-2.642, 2.875, -2.404]), constant_term=-2.443)
    print('equal: {0}'.format(p5 == p6))
    print('parallel but unequal: {0}'.format(p5.is_parallel_to(p6)))

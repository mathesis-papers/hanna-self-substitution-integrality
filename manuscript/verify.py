import json
from fractions import Fraction
from math import comb, factorial, isqrt


DEGREE = 16


def multiply(left, right, degree):
    result = [Fraction(0) for index in range(degree + 1)]
    for left_index, left_value in enumerate(left[:degree + 1]):
        for right_index, right_value in enumerate(right[:degree + 1 - left_index]):
            result[left_index + right_index] += left_value * right_value
    return result


def power(series, exponent, degree):
    result = [Fraction(1)] + [Fraction(0)] * degree
    for index in range(exponent):
        result = multiply(result, series, degree)
    return result


def logarithm_coefficients(series, dilation, weight, degree):
    result = [Fraction(0)] * (degree + 1)
    for outer_index in range(1, degree + 1):
        cutoff = degree - outer_index
        scaled = [value * dilation(outer_index) ** index
                  for index, value in enumerate(series[:cutoff + 1])]
        raised = power(scaled, outer_index, cutoff)
        for inner_index, value in enumerate(raised):
            result[outer_index + inner_index] += Fraction(weight(outer_index), outer_index) * value
    return result


def triangular_solution(dilation, weight, degree):
    series = [Fraction(1)]
    logarithm = [Fraction(0)]
    for current_degree in range(1, degree + 1):
        logarithm.append(logarithm_coefficients(series, dilation, weight, current_degree)[current_degree])
        series.append(sum(index * logarithm[index] * series[current_degree - index]
                          for index in range(1, current_degree + 1)) / current_degree)
    return series


def exponential_by_powers(logarithm, degree):
    result = [Fraction(1)] + [Fraction(0)] * degree
    current_power = result[:]
    for exponent in range(1, degree + 1):
        current_power = multiply(current_power, logarithm, degree)
        result = [value + addition / factorial(exponent)
                  for value, addition in zip(result, current_power)]
    return result


def iterate_solution(dilation, weight, degree):
    series = [Fraction(1)] + [Fraction(0)] * degree
    for iteration in range(degree):
        logarithm = logarithm_coefficients(series, dilation, weight, degree)
        updated = exponential_by_powers(logarithm, degree)
        assert updated[:iteration + 1] == series[:iteration + 1]
        series = updated
    assert exponential_by_powers(logarithm_coefficients(series, dilation, weight, degree), degree) == series
    return series


def sigma(number):
    return sum(divisor for divisor in range(1, number + 1) if number % divisor == 0)


def sparse_solution(base, degree):
    series = [Fraction(1)]
    logarithm = [Fraction(0)]
    for current_degree in range(1, degree + 1):
        logarithm.append(sum(series[current_degree // divisor - 1] * base ** (current_degree - divisor) / divisor
                             for divisor in range(1, current_degree + 1) if current_degree % divisor == 0))
        series.append(sum(index * logarithm[index] * series[current_degree - index]
                          for index in range(1, current_degree + 1)) / current_degree)
    return series


def sparse_iteration(base, degree):
    series = [Fraction(1)] + [Fraction(0)] * degree
    for iteration in range(degree):
        logarithm = [Fraction(0)] * (degree + 1)
        for outer_index in range(1, degree + 1):
            for inner_index in range(degree // outer_index):
                logarithm[outer_index * (inner_index + 1)] += (
                    series[inner_index] * Fraction(base ** (outer_index * inner_index), outer_index))
        series = exponential_by_powers(logarithm, degree)
    return series


def valuation(number, prime):
    exponent = 0
    while number % prime == 0:
        number //= prime
        exponent += 1
    return exponent


def main():
    checks = []
    samples = {
        'A163921': (sigma, lambda index: 1,
                    [1, 1, 2, 7, 41, 385, 5769, 139541, 5551356, 369312953, 41588540350,
                     7987225089655, 2629160183190431, 1487755631073862696, 1450453417949809255147,
                     2439516473122553169216351, 7086426394313598512496200542]),
        'A157675': (lambda index: 2 ** index, lambda index: 1,
                    [1, 1, 3, 19, 237, 5741, 270857, 25099497, 4605241487, 1681614043919,
                     1225216121453227, 1783355695990213771, 5188617952349909215215,
                     30183911091753947571225583, 351131331387346570480797774119]),
        'A159593': (lambda index: 3 ** index, lambda index: 1,
                    [1, 1, 4, 49, 1768, 187474, 58888462, 55210937881, 155033790773008,
                     1305338879106660550, 32966118096763299572020, 2497521410388697783376754490,
                     567627952695201383291867693446222]),
    }
    for sequence_id, (dilation, weight, expected) in samples.items():
        direct = triangular_solution(dilation, weight, DEGREE)
        iterated = iterate_solution(dilation, weight, DEGREE)
        assert direct == iterated
        assert all(value.denominator == 1 for value in direct)
        assert direct[:len(expected)] == expected, sequence_id
        assert all(value >= 0 for value in direct)
        checks.append({'sequence': sequence_id, 'degrees': [0, DEGREE],
                       'source_prefix_degrees': [0, len(expected) - 1],
                       'coefficients': [int(value) for value in direct]})
    shifted_prefixes = {
        2: [1, 1, 2, 18, 476, 38358, 11363548, 15060027956, 92500603618872,
            2483766272252845670, 279689176516909339664044,
            129570236404446129260308225372, 244562582019257683819447274838128648],
        3: [1, 1, 3, 72, 10899, 13775400, 185126087574, 25804839804745266,
            34576186149297464460939, 431016348991655088520874274528,
            49345807665910247245278532420622266770],
    }
    for base in [-2, -1, 0, 1, 2, 3]:
        dilation = lambda index: base ** index
        weight = lambda index: base ** (index * index)
        direct = triangular_solution(dilation, weight, DEGREE - 1)
        assert direct == iterate_solution(dilation, weight, DEGREE - 1)
        shifted = [Fraction(1)] + direct
        assert all(value.denominator == 1 for value in shifted)
        expected = shifted_prefixes.get(base)
        if expected:
            assert shifted[:len(expected)] == expected
            assert all(value >= 0 for value in shifted)
        checks.append({'family': 'shifted', 'base': base, 'degrees': [0, DEGREE],
                       'coefficients': [int(value) for value in shifted]})
    catalan = triangular_solution(lambda index: 1, lambda index: 1, DEGREE)
    assert catalan == [Fraction(comb(2 * index, index), index + 1) for index in range(DEGREE + 1)]
    weak = triangular_solution(lambda index: 2 if index % 4 == 0 else 0, lambda index: 1, DEGREE)
    assert weak == iterate_solution(lambda index: 2 if index % 4 == 0 else 0, lambda index: 1, DEGREE)
    assert all(value.denominator == 1 for value in weak)
    sparse = sparse_solution(2, DEGREE)
    assert sparse == sparse_iteration(2, DEGREE)
    expected_sparse = [1, 1, 3, 15, 139, 2387, 79115, 5148411, 664332843, 170744863371,
                       87593505706987, 89783692196468907, 183966962290186844267,
                       753712824966410639243755, 6175169543791440589003293035,
                       101180154484297968338398947674219]
    assert sparse[:len(expected_sparse)] == expected_sparse
    assert all(value.denominator == 1 for value in sparse)
    assert all(value >= 0 for value in sparse)
    checks.append({'sequence': 'A179471', 'degrees': [0, DEGREE],
                   'source_prefix_degrees': [0, len(expected_sparse) - 1],
                   'coefficients': [int(value) for value in sparse]})
    for base in [-2, -1, 0, 1, 3]:
        direct = sparse_solution(base, DEGREE)
        assert direct == sparse_iteration(base, DEGREE)
        assert all(value.denominator == 1 for value in direct)
        if base >= 0:
            assert all(value >= 0 for value in direct)
        checks.append({'family': 'sparse', 'base': base, 'degrees': [0, DEGREE],
                       'coefficients': [int(value) for value in direct]})
    euler_exponents = [-2, 1, 0, 3, -1, 2]
    weighted_inputs = {
        'divisor-sum weights': sigma,
        'signed geometric weights': lambda index: (-2) ** index,
        'signed Euler weights': lambda index: sum(
            divisor * euler_exponents[divisor - 1]
            for divisor in range(1, len(euler_exponents) + 1) if index % divisor == 0),
    }
    for name, weight in weighted_inputs.items():
        direct = triangular_solution(sigma, weight, DEGREE)
        assert direct == iterate_solution(sigma, weight, DEGREE)
        assert all(value.denominator == 1 for value in direct)
        checks.append({'family': 'two integral inputs', 'weight': name,
                       'degrees': [0, DEGREE],
                       'coefficients': [int(value) for value in direct]})
    for prime in [2, 3]:
        dilation = lambda index: prime if index % (prime * prime) == 0 else 0
        direct = triangular_solution(dilation, lambda index: 1, DEGREE)
        assert direct == iterate_solution(dilation, lambda index: 1, DEGREE)
        assert all(value.denominator == 1 for value in direct)
        logarithm = [Fraction(0)] + [Fraction(dilation(index), index)
                                     for index in range(1, DEGREE + 1)]
        input_series = exponential_by_powers(logarithm, DEGREE)
        assert input_series[prime * prime] == Fraction(1, prime)
        checks.append({'family': 'converse obstruction not asserted in source',
                       'prime': prime, 'degrees': [0, DEGREE],
                       'fixed_point_coefficients': [int(value) for value in direct],
                       'input_coefficient_degree': prime * prime,
                       'input_coefficient': str(input_series[prime * prime])})
    rational_control = triangular_solution(lambda index: Fraction(1, 2), lambda index: 1, 4)
    assert rational_control == iterate_solution(lambda index: Fraction(1, 2), lambda index: 1, 4)
    assert rational_control[2] == Fraction(3, 2)
    assert sparse[0] == 1
    primes = [candidate for candidate in range(2, 32)
              if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1))]
    congruence_count = 0
    for prime in primes:
        for number in range(1, 41):
            modulus = prime ** (valuation(number, prime) + 1)
            for base in range(-3, 5):
                assert (pow(base, (prime * number) ** 2, modulus) - pow(base, number ** 2, modulus)) % modulus == 0
                assert (pow(base, prime * number, modulus) - pow(base, number, modulus)) % modulus == 0
                congruence_count += 2
    result = {'status': 'PASS', 'arithmetic': 'exact rational integers; no rounding',
              'scope': 'finite coefficient identities and congruences', 'series': checks,
              'additional_checks': ['Catalan formula, degrees 0..16',
                                    'weaker hypothesis example, degrees 0..16',
                                    'all six source prefixes have nonnegative coefficients',
                                    'unrestricted rational dilation control: degree 2 coefficient 3/2',
                                    'A179471 printed logarithm has constant coefficient 1; corrected logarithm has 0'],
              'congruences': {'prime_bound_inclusive': 31, 'index_bound_inclusive': 40,
                              'bases_inclusive': [-3, 4], 'count': congruence_count}}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

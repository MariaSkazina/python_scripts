import unittest

from fastq_barbershop import calculate_gc_content, filter,\
    leading, sliding_window, trailing, symbols2values, crop, headcrop


class FilterTest(unittest.TestCase):
    def test_calculate_gc_content_empty_read(self):
        self.assertEqual(calculate_gc_content(''), 0)

    def test_calculate_gc_content_only_G(self):
        self.assertEqual(calculate_gc_content('ATTG'), 25)

    def test_calculate_gc_content_only_C(self):
        self.assertEqual(calculate_gc_content('ATCT'), 25)

    def test_calculate_gc_content_normal(self):
        self.assertEqual(calculate_gc_content('ATGGGCCTAA'), 50)

    def test_filter_less_min_length(self):
        self.assertEqual(filter('ATGC', 10, 0, 100), 'bad')

    def test_filter_normal_length(self):
        self.assertEqual(filter('ATGC', 4, 0, 100), 'good')

    def test_filter_less_min_gc(self):
        self.assertEqual(filter('ATGAACTTTAA', 0, 25, 75), 'bad')

    def test_filter_normal_gc(self):
        self.assertEqual(filter('ATGGCGACT', 0, 25, 75), 'good')

    def test_filter_more_max_gc(self):
        self.assertEqual(filter('GCGACGTGC', 0, 25, 75), 'bad')


class TrimmomaticTest(unittest.TestCase):

    def test_symbol2values_check_quality(self):
        self.assertEqual(symbols2values('!!!'), [0, 0, 0])

    def test_leading_not_cut(self):
        self.assertEqual(leading('ACCCGT', 'AAABBA', 10), ('ACCCGT', 'AAABBA'))

    def test_leading_cut(self):
        self.assertEqual(leading('AAA', '!!!', 42), ('', ''))

    def test_leading_cut_not_all(self):
        self.assertEqual(leading('ACG', '!AA', 30), ('CG', 'AA'))

    def test_leading_not_cut_from_other_side(self):
        self.assertEqual(leading('ACGACGACG', 'AAAAAAA!!', 30), ('ACGACGACG', 'AAAAAAA!!'))

    def test_trailing_not_cut(self):
        self.assertEqual(trailing('ACCCGT', 'AAABBA', 10), ('ACCCGT', 'AAABBA'))

    def test_trailing_cut(self):
        self.assertEqual(trailing('AAA', '!!!', 10), ('', ''))

    def test_trailing_cut_not_all(self):
        self.assertEqual(trailing('ACG', 'AA!', 30), ('AC', 'AA'))

    def test_trailing_not_cut_from_other_side(self):
        self.assertEqual(trailing('ACGACGACG', '!!AAAAAAA', 30), ('ACGACGACG', '!!AAAAAAA'))

    def test_sliding_window_not_cut(self):
        self.assertEqual(sliding_window('CAGTCAGTC', 'AAAAAAAAA', 32, 2), ('CAGTCAGTC', 'AAAAAAAAA'))

    def test_sliding_window_cut_all(self):
        self.assertEqual(sliding_window('CAGTCAGTC', 'AAA3567@?', 32, 4), ('', ''))

    def test_sliding_window_cut_dramatic_quality_drop(self):
        self.assertEqual(sliding_window('CAGTCAGTC', 'AAAAAA!!!', 32, 3), ('CAGT', 'AAAA'))

    def test_sliding_window_cut_not_so_dramatic_quality_drop(self):
        self.assertEqual(sliding_window('CAGTCAGTC', 'AAAAAAA@?', 32, 2), ('CAGTCA', 'AAAAAA'))

    def test_headcrop_crop(self):
        self.assertEqual(headcrop('AGTCAGTCAGTC', 'ABCABCABCABC', 3), ('CAGTCAGTC', 'ABCABCABC'))

    def test_headcrop_notcrop(self):
        self.assertEqual(headcrop('AGTCAGTCAGTC', 'ABCABCABCABC', 0), ('AGTCAGTCAGTC', 'ABCABCABCABC'))

    def test_headcrop_length_more_than_read(self):
        self.assertRaises(ValueError, lambda: headcrop('AGTCAGTCAGTC', 'AABCDABCDEFA', 20))

    def test_crop_crop(self):
        self.assertEqual(crop('AGTCAGTCAGTG', 'ABCABCABCABC',  3), ('AGTCAGTCA', 'ABCABCABC'))

    def test_crop_notcrop(self):
        self.assertEqual(crop('AGTCAGTCAGTC', 'ABCABCABCABC', 0), ('AGTCAGTCAGTC', 'ABCABCABCABC'))

    def test_crop_length_more_than_read(self):
        self.assertRaises(ValueError, lambda: crop('AGTCAGTCAGTC', 'ABCABCABCABC', 20))


if __name__ == '__main__':
    unittest.main()

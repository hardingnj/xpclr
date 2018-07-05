from xpclr import util


def test_load_vcf():

    path = 'fixture/small.vcf.gz'

    samplesA = "AB0085-C,AB0087-C,AB0088-C,AB0089-C,AB0090-C"
    samplesB = "AN0011-C,AN0012-C,AN0014-C,AN0016-C,AN0017-C,AN0018-C"

    g1, g2, pos, gdist = util.load_vcf_format_data(path, chrom="3L", s1=samplesA, s2=samplesB)

    assert g1.shape[0] == g2.shape[0] == pos.shape[0]
    assert g1.shape[1] == 5
    assert g2.shape[1] == 6

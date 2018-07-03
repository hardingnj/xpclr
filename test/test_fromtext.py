from xpclr import util

def test_loadtext():

    g1, g2, positions, genetic_dist = util.load_text_format_data(
        "fixture/mapfile.snp",
        "fixture/genotype1.geno",
        "fixture/genotype2.geno")

    assert g1.shape[0] == 200
    assert g2.shape[0] == 200

    assert positions.dtype == "int"
    assert genetic_dist.dtype == "float"



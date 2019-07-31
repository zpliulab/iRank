import cPickle
import pkg_resources
#con1 =open("pheno_cooccur.pickle","rb")
dir_data="data_phenorank"
con = open(pkg_resources.resource_filename("phenorank", dir_data + "/pheno_cooccur.pickle"), "rb")
gc_h = cPickle.load(con)
gc_h

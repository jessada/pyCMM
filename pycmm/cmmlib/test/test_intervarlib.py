import unittest
from pycmm.template import SafeTester
from pycmm.cmmlib.intervarlib import parse_intervar_class
from pycmm.cmmlib.intervarlib import parse_intervar_evidence
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_PATHOGENIC
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_PATHOGENIC


class TestDNARegion(SafeTester):

    def __init__(self, methodName):
        super(TestDNARegion, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def test_parse_intervar_0(self):
        """ test blank intervar """

        self.init_test(self.current_func_name)
        raw_intervar = ""
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         "",
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def test_parse_intervar_1(self):
        """ test benign intervar """

        self.init_test(self.current_func_name)
        raw_intervar = "InterVar:Benign;PVS1=0;PS=[0;0;0;0;0];PM=[0;0;0;0;0;0;0];PP=[0;0;0;0;0;0];BA1=1;BS=[1;1;0;0;0];BP=[0;0;0;1;0;1;0;0]"
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         INTERVAR_CLASS_BENIGN,
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "BA1, BS1, BS2, BP4, BP6",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def test_parse_intervar_2(self):
        """ test likely-benign intervar """

        self.init_test(self.current_func_name)
        raw_intervar = "InterVar:Likelybenign;PVS1=0;PS=[0;0;0;0;0];PM=[0;0;0;0;0;0;0];PP=[0;0;0;0;0;0];BA1=0;BS=[0;1;0;0;0];BP=[1;0;0;0;0;1;0;0]"
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         INTERVAR_CLASS_LIKELY_BENIGN,
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "BS2, BP1, BP6",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def test_parse_intervar_3(self):
        """ test uncertain-significance intervar """

        self.init_test(self.current_func_name)
        raw_intervar = "InterVar:UncertainSignificance;PVS1=0;PS=[0;0;0;0;0];PM=[0;1;0;0;0;0;0];PP=[0;0;0;0;0;0];BA1=0;BS=[0;1;0;0;0];BP=[0;0;0;0;0;0;0;0]"
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE,
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "PM2, BS2",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def test_parse_intervar_4(self):
        """ test likely-pathogenic intervar """

        self.init_test(self.current_func_name)
        raw_intervar = "InterVar:Likelypathogenic;PVS1=1;PS=[0;0;0;0;0];PM=[0;1;0;0;0;0;0];PP=[0;0;0;0;0;0];BA1=0;BS=[0;0;0;0;0];BP=[0;0;0;0;0;0;0;0]"
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         INTERVAR_CLASS_LIKELY_PATHOGENIC,
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "PVS1, PM2",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def test_parse_intervar_5(self):
        """ test pathogenic intervar """

        self.init_test(self.current_func_name)
        raw_intervar = "InterVar:Pathogenic;PVS1=1;PS=[0;0;0;0;0];PM=[0;1;0;0;0;0;0];PP=[0;0;0;0;1;0];BA1=0;BS=[0;0;0;0;0];BP=[0;0;0;0;0;0;0;0]"
        intervar_class = parse_intervar_class(raw_intervar)
        self.assertEqual(intervar_class,
                         INTERVAR_CLASS_PATHOGENIC,
                         "Cannot correctly parse intervar: " + intervar_class)
        intervar_evidence = parse_intervar_evidence(raw_intervar)
        self.assertEqual(intervar_evidence,
                         "PVS1, PM2, PP5",
                         "Cannot correctly parse intervar: " + str(intervar_evidence))

    def tearDown(self):
        self.remove_working_dir()

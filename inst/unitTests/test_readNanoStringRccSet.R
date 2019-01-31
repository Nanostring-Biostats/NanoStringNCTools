test_readNanoStringRccSet_onearg <- function() {
  datadir <- system.file("extdata", "3D_Bio_Example_Data",
                         package = "NanoStringNCTools")
  rcc <-
    readNanoStringRccSet(dir(datadir, pattern = "SKMEL.*\\.RCC$",
                             full.names = TRUE))
  checkTrue(validObject(rcc))
  checkIdentical(c(Features = 397L, Samples = 12L), dim(rcc))
  checkIdentical(data.frame(BarcodeClass = c("Endogenous", "SNV_REF"),
                            GeneName =
                              c("TP53",
                                "PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141"),
                            Accession = c("NM_000546.2", "nRef_00032.1"),
                            IsControl = c(FALSE, FALSE),
                            ControlConc = c(NA_real_, NA_real_),
                            row.names =
                              c("Endogenous_TP53_NM_000546.2",
                                "SNV_REF_PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141_nRef_00032.1"),
                            stringsAsFactors = FALSE),
                 fData(rcc)[c(1L, nrow(rcc)),])
  checkIdentical(data.frame(FileVersion = numeric_version(c("1.7", "1.7")),
                            SoftwareVersion = numeric_version(c("4.0.0.3", "4.0.0.3")),
                            SystemType = c("Gen2", "Gen2"),
                            SampleID = c("SKMEL2-DMSO-8hr", "SKMEL28-V-8hr"),
                            SampleOwner = c("", ""),
                            SampleComments = c("DNA-RNA-Protein", "DNA-RNA-Protein"),
                            SampleDate = as.Date(c("2017-01-13", "2017-01-25")),
                            SystemAPF = c("n6_vDV1", "n6_vDV1"),
                            AssayType = c(NA_character_, NA_character_),
                            LaneID = c(4L, 10L),
                            FovCount = c(280L, 280L),
                            FovCounted = c(268L, 279L),
                            ScannerID = c("1207C0049", "1207C0049"),
                            StagePosition = c(3L, 3L),
                            BindingDensity = c(0.8, 0.64),
                            CartridgeID = c("AGBT-3DBio-C1-SKMEL2",
                                            "AGBT-3DBio-2-Repeat-C3-SKMEL28"),
                            CartridgeBarcode = c("", ""),
                            row.names = c("SKMEL2-DMSO-8h-R1_04.RCC",
                                          "SKMEL28-VEM-8h-R3_10.RCC"),
                            stringsAsFactors = FALSE),
                 sData(rcc)[c(1L, ncol(rcc)),])
}

test_readNanoStringRccSet_rlf <- function() {
  datadir <- system.file("extdata", "3D_Bio_Example_Data",
                         package = "NanoStringNCTools")
  rcc <-
    readNanoStringRccSet(dir(datadir, pattern = "SKMEL.*\\.RCC$",
                             full.names = TRUE),
                         file.path(datadir, "3D_SolidTumor_Sig.rlf"))
  checkTrue(validObject(rcc))
  checkIdentical(c(Features = 397L, Samples = 12L), dim(rcc))
  checkIdentical(data.frame(BarcodeClass = c("Endogenous", "SNV_REF"),
                            GeneName =
                              c("TP53",
                                "PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141"),
                            Accession = c("NM_000546.2", "nRef_00032.1"),
                            IsControl = c(FALSE, FALSE),
                            ControlConc = c(NA_real_, NA_real_),
                            Barcode = c("BGBGYR", "YRGRYR"),
                            ProbeID = c("NM_000546.2:1330", "nRef_00032:1"),
                            Species = c("Hs", "Hs"),
                            BarcodeComments = c("ENDOGENOUS", "SNV_REF"),
                            row.names =
                              c("Endogenous_TP53_NM_000546.2",
                                "SNV_REF_PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141_nRef_00032.1"),
                            stringsAsFactors = FALSE),
                 fData(rcc)[c(1L, nrow(rcc)),
                            c("BarcodeClass", "GeneName", "Accession",
                              "IsControl", "ControlConc",
                              "Barcode", "ProbeID", "Species",
                              "BarcodeComments")])
  checkIdentical(data.frame(FileVersion = numeric_version(c("1.7", "1.7")),
                            SoftwareVersion = numeric_version(c("4.0.0.3", "4.0.0.3")),
                            SystemType = c("Gen2", "Gen2"),
                            SampleID = c("SKMEL2-DMSO-8hr", "SKMEL28-V-8hr"),
                            SampleOwner = c("", ""),
                            SampleComments = c("DNA-RNA-Protein", "DNA-RNA-Protein"),
                            SampleDate = as.Date(c("2017-01-13", "2017-01-25")),
                            SystemAPF = c("n6_vDV1", "n6_vDV1"),
                            AssayType = c(NA_character_, NA_character_),
                            LaneID = c(4L, 10L),
                            FovCount = c(280L, 280L),
                            FovCounted = c(268L, 279L),
                            ScannerID = c("1207C0049", "1207C0049"),
                            StagePosition = c(3L, 3L),
                            BindingDensity = c(0.8, 0.64),
                            CartridgeID = c("AGBT-3DBio-C1-SKMEL2",
                                            "AGBT-3DBio-2-Repeat-C3-SKMEL28"),
                            CartridgeBarcode = c("", ""),
                            row.names = c("SKMEL2-DMSO-8h-R1_04.RCC",
                                          "SKMEL28-VEM-8h-R3_10.RCC"),
                            stringsAsFactors = FALSE),
                 sData(rcc)[c(1L, ncol(rcc)),])
}

test_readNanoStringRccSet_rlf_pheno <- function() {
  datadir <- system.file("extdata", "3D_Bio_Example_Data",
                         package = "NanoStringNCTools")
  rcc <-
    readNanoStringRccSet(dir(datadir, pattern = "SKMEL.*\\.RCC$",
                             full.names = TRUE),
                         file.path(datadir, "3D_SolidTumor_Sig.rlf"),
                         file.path(datadir, "3D_SolidTumor_PhenoData.csv"))
  checkTrue(validObject(rcc))
  checkIdentical(c(Features = 397L, Samples = 12L), dim(rcc))
  checkIdentical(data.frame(BarcodeClass = c("Endogenous", "SNV_REF"),
                            GeneName =
                              c("TP53",
                                "PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141"),
                            Accession = c("NM_000546.2", "nRef_00032.1"),
                            Barcode = c("BGBGYR", "YRGRYR"),
                            ProbeID = c("NM_000546.2:1330", "nRef_00032:1"),
                            Species = c("Hs", "Hs"),
                            BarcodeComments = c("ENDOGENOUS", "SNV_REF"),
                            IsControl = c(FALSE, FALSE),
                            row.names =
                              c("Endogenous_TP53_NM_000546.2",
                                "SNV_REF_PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141_nRef_00032.1"),
                            stringsAsFactors = FALSE),
                 fData(rcc)[c(1L, nrow(rcc)),
                            c("BarcodeClass", "GeneName", "Accession",
                              "Barcode", "ProbeID", "Species",
                              "BarcodeComments", "IsControl")])
  checkIdentical(data.frame(Treatment = c("DMSO", "VEM"),
                            BRAFGenotype = c("wt/wt", "mut/mut"),
                            FileVersion = numeric_version(c("1.7", "1.7")),
                            SoftwareVersion = numeric_version(c("4.0.0.3", "4.0.0.3")),
                            SystemType = c("Gen2", "Gen2"),
                            SampleID = c("SKMEL2-DMSO-8hr", "SKMEL28-V-8hr"),
                            SampleOwner = c("", ""),
                            SampleComments = c("DNA-RNA-Protein", "DNA-RNA-Protein"),
                            SampleDate = as.Date(c("2017-01-13", "2017-01-25")),
                            SystemAPF = c("n6_vDV1", "n6_vDV1"),
                            AssayType = c(NA_character_, NA_character_),
                            LaneID = c(4L, 10L),
                            FovCount = c(280L, 280L),
                            FovCounted = c(268L, 279L),
                            ScannerID = c("1207C0049", "1207C0049"),
                            StagePosition = c(3L, 3L),
                            BindingDensity = c(0.8, 0.64),
                            CartridgeID = c("AGBT-3DBio-C1-SKMEL2",
                                            "AGBT-3DBio-2-Repeat-C3-SKMEL28"),
                            CartridgeBarcode = c("", ""),
                            row.names = c("SKMEL2-DMSO-8h-R1_04.RCC",
                                          "SKMEL28-VEM-8h-R3_10.RCC"),
                            stringsAsFactors = FALSE),
                 sData(rcc)[c(1L, ncol(rcc)),])
}

test_readNanoStringRccSet_rlf_pheno <- function() {
  datadir <- system.file("extdata", "3D_Bio_Example_Data",
                         package = "NanoStringNCTools")
  rcc <-
    readNanoStringRccSet(dir(datadir, pattern = "SKMEL.*\\.RCC$",
                             full.names = TRUE),
                         file.path(datadir, "3D_SolidTumor_Sig.rlf"),
                         file.path(datadir, "3D_SolidTumor_PhenoData.csv"),
                         phenoDataColPrefix = "PHENO_")
  checkTrue(validObject(rcc))
  checkIdentical(c(Features = 397L, Samples = 12L), dim(rcc))
  checkIdentical(data.frame(BarcodeClass = c("Endogenous", "SNV_REF"),
                            GeneName =
                              c("TP53",
                                "PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141"),
                            Accession = c("NM_000546.2", "nRef_00032.1"),
                            Barcode = c("BGBGYR", "YRGRYR"),
                            ProbeID = c("NM_000546.2:1330", "nRef_00032:1"),
                            Species = c("Hs", "Hs"),
                            BarcodeComments = c("ENDOGENOUS", "SNV_REF"),
                            IsControl = c(FALSE, FALSE),
                            row.names =
                              c("Endogenous_TP53_NM_000546.2",
                                "SNV_REF_PIK3CA Ref (exon 10)|hg19|+|chr3:178936060-178936141_nRef_00032.1"),
                            stringsAsFactors = FALSE),
                 fData(rcc)[c(1L, nrow(rcc)),
                            c("BarcodeClass", "GeneName", "Accession",
                              "Barcode", "ProbeID", "Species",
                              "BarcodeComments", "IsControl")])
  checkIdentical(data.frame(PHENO_Treatment = c("DMSO", "VEM"),
                            PHENO_BRAFGenotype = c("wt/wt", "mut/mut"),
                            FileVersion = numeric_version(c("1.7", "1.7")),
                            SoftwareVersion = numeric_version(c("4.0.0.3", "4.0.0.3")),
                            SystemType = c("Gen2", "Gen2"),
                            SampleID = c("SKMEL2-DMSO-8hr", "SKMEL28-V-8hr"),
                            SampleOwner = c("", ""),
                            SampleComments = c("DNA-RNA-Protein", "DNA-RNA-Protein"),
                            SampleDate = as.Date(c("2017-01-13", "2017-01-25")),
                            SystemAPF = c("n6_vDV1", "n6_vDV1"),
                            AssayType = c(NA_character_, NA_character_),
                            LaneID = c(4L, 10L),
                            FovCount = c(280L, 280L),
                            FovCounted = c(268L, 279L),
                            ScannerID = c("1207C0049", "1207C0049"),
                            StagePosition = c(3L, 3L),
                            BindingDensity = c(0.8, 0.64),
                            CartridgeID = c("AGBT-3DBio-C1-SKMEL2",
                                            "AGBT-3DBio-2-Repeat-C3-SKMEL28"),
                            CartridgeBarcode = c("", ""),
                            row.names = c("SKMEL2-DMSO-8h-R1_04.RCC",
                                          "SKMEL28-VEM-8h-R3_10.RCC"),
                            stringsAsFactors = FALSE),
                 sData(rcc)[c(1L, ncol(rcc)),])
}

rm(list = ls())
gc()
setwd("Path to data")

###导入模块，输入gene id
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)

# #FL
# gene_id<-c('ANKRD26P1', 'C1orf109', 'ZNF544', 'CDKN2C', 'SHTN1', 'DLGAP1', 'HMBOX1', 'ARFGAP1', 'GTSE1', 'CDCA7L', 'UBTD2', 'ANP32B', 'DMTN', 'COG8', 'MCM2', 'PCGF6', 'CCDC102B', 'KRTAP10-11', 'KRTAP10-3', 'KRTAP10-1', 'KRTAP10-9', 'KRTAP10-2', 'KRTAP10-4', 'KRTAP10-12', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-10', 'KRTAP10-5', 'FAAP20', 'CIB2', 'ZNF572', 'C8orf34', 'NUFIP2', 'ACTN3', 'ACTN4', 'KLHL38', 'SPRED2', 'TRIM10', 'CBWD1', 'CBWD3', 'CBWD5', 'CBWD6', 'CBWD2', 'UBXN2B', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G3', 'OR7G2', 'OR7D2', 'OR7G1', 'POMGNT1', 'HECTD3', 'PRKAG1', 'GNL3', 'MALSU1', 'PPP2CA', 'MYL3', 'COL20A1', 'RB1', 'TCEANC', 'CCNH', 'KYNU', 'SPAG8', 'BEND2', 'MID2', 'MRFAP1L1', 'MRFAP1', 'MCMBP', 'CPNE6', 'ARNTL', 'CNPY2', 'FXR2', 'TXNDC9', 'CCNC', 'EMSY', 'PPP6R2', 'BLOC1S5-TXNDC5', 'TXNDC5')
# #tox8
# gene_id<-c('RHOA', 'TADA2B', 'UBTD2', 'EIF2B2', 'CNST', 'CKB', 'CKM', 'UBE2R2', 'MDM2', 'KAZN', 'EMSY', 'NEIL2', 'CCNC', 'FXR2')
# #tox7
# gene_id<-c('LYSMD1', 'ZMIZ2', 'GTSE1', 'UBTD2', 'EMSY', 'G2E3', 'TADA2B', 'HEPACAM', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G2', 'OR7G3', 'OR7E37P', 'OR7E5P', 'OR7D2', 'OR7E12P', 'OR7E14P', 'OR7E47P', 'OR7E2P', 'OR7G1', 'OR1M1', 'SFN', 'EIF2B2', 'AP4M1', 'CKB', 'CKM', 'PLCD4', 'VCL', 'EML1', 'ACD', 'FANCL', 'ESRRG', 'TRIM27', 'CNST', 'NTM', 'OPCML', 'MORN3', 'CNPY2', 'FXR2', 'KAZN', 'CCNC', 'NEIL2', 'BLOC1S5-TXNDC5', 'TXNDC5', 'MDM2')
# #tox6
# gene_id<-c('PRKAA2', 'UBTD2', 'G2E3', 'SNAPC5', 'TADA2B', 'NUFIP2', 'ZMIZ2', 'XRCC4', 'MLLT11', 'SFN', 'RHOA', 'ZNF23', 'CKB', 'CKM', 'CNPY2', 'FXR2', 'CCNC', 'EMSY', 'BLOC1S5-TXNDC5', 'TXNDC5')
# #tox5
# gene_id<-c('KANK2', 'HMBOX1', 'GTSE1', 'RBP3', 'TPM3', 'TPM3P9', 'TPM4', 'DUS1L', 'DMTN', 'DTX2', 'QARS1', 'SEC23B', 'SNRNP25', 'G2E3', 'ING3', 'SNAPC5', 'PCNP', 'COG4', 'NTM', 'OPCML', 'KRTAP10-3', 'KRTAP10-9', 'KRTAP10-1', 'KRTAP10-11', 'KRTAP10-2', 'KRTAP10-4', 'KRTAP10-10', 'KRTAP10-12', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-5', 'ACD', 'CCDC84', 'CPSF1', 'PGGHG', 'CCL16', 'BAGE2', 'BAGE3', 'BAGE5', 'BAGE4', 'BAGE', 'GCM2', 'PDE4C', 'PDE4A', 'FAM74A4', 'FAM74A6', 'FAM74A3', 'FAM74A7', 'FAM74A1', 'RBBP7', 'NELL1', 'TADA2B', 'PHF12', 'UTP14C', 'ALG11', 'UTP14A', 'STAC2', 'CD22', 'DGKG', 'CDX1', 'SNAPC3', 'SMCO1', 'CLIP1', 'MON2', 'NUFIP2', 'CERT1', 'ACTN3', 'ACTN4', 'SMARCB1', 'UBXN2B', 'TMEM129', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G3', 'OR7D2', 'OR7G2', 'OR7E14P', 'OR7G1', 'XRCC4', 'KIAA0825', 'FCAR', 'HECTD3', 'PRKAG1', 'MALSU1', 'PPP2CA', 'MYL3', 'HDAC7', 'UEVLD', 'PPP1R8', 'XRN2', 'WASHC3', 'C1orf116', 'COPS8', 'SLA', 'RHOA', 'CNPY2', 'SCNM1', 'TNFAIP8L2-SCNM1', 'RRS1', 'EMSY', 'PRKAG2', 'C1QTNF1', 'CD300LG', 'ZNF23', 'ZNF420', 'ZNF619', 'ZNF487', 'ZNF700', 'ZNF419', 'ZNF773', 'ZNF772', 'ZNF304', 'ZIK1', 'ZNF551', 'ZNF814', 'ZNF776', 'ZNF587', 'ZNF417', 'ZNF418', 'CAGE1', 'NAV1', 'SIGLEC6', 'SIGLEC17P', 'CD33', 'KDM8', 'INPP5K', 'ZBTB16', 'ZNF202', 'TRIP10', 'RIMBP2', 'TUBG1', 'TUBG2', 'ARIH2', 'USB1', 'FAM89A', 'CEP19', 'RFPL3', 'RFPL2', 'RFPL1', 'RFPL3S', 'RFPL1S', 'PSMD2', 'ZSWIM2', 'PACSIN2', 'MLX', 'PRKACB', 'MAPRE3', 'C1QTNF6', 'STARD7', 'OTUB1', 'CMBL', 'XPO5', 'PBX3', 'SNF8', 'NIBAN1', 'NIFK', 'CCNH', 'CIAPIN1', 'GLRX3', 'KYNU', 'EME1', 'ASH2L', 'SLC22A17', 'PDE4D', 'MINDY1', 'INTS12', 'GABRG1', 'NR1H2', 'SNPH', 'ZC3HC1', 'NELFA', 'SPAG8', 'CPNE4', 'SMAD1', 'SMAD1-AS1', 'GATA2', 'PNKP', 'TSTD2', 'NOP58', 'GRN', 'CAPN11', 'NELFCD', 'LILRB1', 'LILRA3', 'LILRA1', 'LILRA2', 'LILRB2', 'LILRB4', 'LILRB3', 'LILRB5', 'EFHC1', 'VPS33A', 'PLCD4', 'CFAP70', 'MTMR4', 'PACSIN1', 'TRIM32', 'GP1BA', 'TGM2', 'EML1', 'RNH1', 'MRFAP1L1', 'MRFAP1', 'MCMBP', 'CPNE6', 'DNAI2', 'ESRRG', 'EYA3', 'DRC3', 'ARNTL', 'KIAA1586', 'PARD6B', 'SCAF1', 'MORN3', 'ZNF410', 'ASPH', 'NMRK1', 'CTDSP2', 'CTDSPL', 'FXR2', 'DCUN1D4', 'CCDC103', 'FAM221A', 'MKRN3', 'BLOC1S5-TXNDC5', 'TXNDC5', 'CRYBA1', 'CIB3', 'STRIP2', 'RNF181', 'MID2', 'CDH11', 'HGS', 'CA8', 'PRMT7')
# #tox4
# gene_id<-c('C1orf109', 'RWDD4', 'POPDC2', 'THADA', 'ARPC3', 'CARD9', 'SHTN1', 'KLHL32', 'CIR1', 'VPS39', 'KANK2', 'ZFYVE26', 'HMBOX1', 'ARFGAP1', 'CDCA7L', 'TPM3', 'TPM3P9', 'TPM4', 'POMP', 'CYREN', 'DMTN', 'TAF1B', 'LMBR1L', 'C4orf33', 'RABL6', 'DTX2', 'GYS1', 'CHCHD2', 'CNPY2', 'G2E3', 'ING3', 'EIF4EBP1', 'SNAPC5', 'NUP62CL', 'PPP2CA', 'ZNF138', 'ZNF92', 'ZNF738', 'ZNF254', 'ZNF506', 'ZNF99', 'ZNF430', 'ZNF676', 'ZNF493', 'ZNF626', 'ZNF708', 'ZNF681', 'ZNF729', 'ZNF100', 'ZNF675', 'ZNF91', 'ZNF107', 'ZNF429', 'ZNF678', 'ZNF90', 'ZNF117', 'ERV3-1-ZNF117', 'ZNF730', 'ZNF695', 'ZNF718', 'ZNF66', 'ZNF732', 'KRTAP10-11', 'KRTAP10-9', 'KRTAP10-1', 'KRTAP10-2', 'KRTAP10-5', 'KRTAP10-4', 'KRTAP10-12', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-3', 'KRTAP10-10', 'KRTAP10-8', 'LCN8', 'CPSF1', 'ZNF572', 'PGGHG', 'PDE4C', 'PDE4A', 'FAM74A4', 'FAM74A6', 'FAM74A3', 'FAM74A7', 'FAM74A1', 'CKB', 'CKM', 'ZNF844', 'ZNF20', 'ZNF625-ZNF20', 'ZNF625', 'ZNF440', 'ZNF700', 'ZNF439', 'ZNF878', 'ZNF799', 'ZNF627', 'ZNF136', 'ZNF443', 'ZNF563', 'ZNF669', 'ZNF101', 'ZNF709', 'GPR68', 'RTN1', 'NELL1', 'RIPPLY1', 'PMS1', 'UTP14C', 'ALG11', 'UTP14A', 'ZBTB49', 'STAC2', 'DGKG', 'CDX1', 'SNAPC3', 'CCDC198', 'UBXN2B', 'OIT3', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G2', 'OR7D2', 'OR7G3', 'OR7E37P', 'OR7E5P', 'OR7E14P', 'OR7E47P', 'OR7E2P', 'OR7G1', 'SFI1', 'KIAA0825', 'CDC34', 'THAP10', 'UBTD1', 'SMYD3', 'HOXA1', 'ING5', 'CUTC', 'HECTD3', 'PRKAG1', 'COL4A3BP', 'SLA', 'COQ3', 'YJU2', 'TRAPPC2L', 'EML1', 'RASA1', 'MCM2', 'MCFD2', 'USHBP1', 'NMRK1', 'COPB1', 'COL20A1', 'ATG12', 'TTC14', 'RBX1', 'ZNF593', 'CARF', 'LYNX1-SLURP2', 'C11orf49', 'CFAP298', 'C21orf59-TCP10L', 'RRS1', 'CINP', 'CPSF4', 'GFOD2', 'UBA6', 'WDR61', 'ZNF23', 'ZNF420', 'ZNF619', 'ZNF487', 'CAGE1', 'SIGLEC6', 'SIGLEC17P', 'CD33', 'ZBTB16', 'ZNF202', 'RIMBP2', 'BANP', 'SNX29P2', 'GPANK1', 'CCNH', 'TBC1D23', 'BABAM1', 'GLRX3', 'KYNU', 'NEK6', 'NAGK', 'NTN4', 'ZNF830', 'ASH2L', 'CCNC', 'PDE4D', 'SNPH', 'MCRS1', 'SMAD1', 'SMAD1-AS1', 'PNKP', 'MBD4', 'TSTD2', 'MID2', 'LCK', 'PACSIN2', 'PACSIN1', 'FAM129A', 'STAC', 'RNH1', 'AHNAK', 'MRFAP1L1', 'MRFAP1', 'CPNE6', 'ACD', 'ASPH', 'RHOJ', 'ATG9A', 'RBBP5', 'KAZN', 'PPP6R2', 'BLOC1S5-TXNDC5', 'TXNDC5', 'KYAT1', 'CCDC103', 'FAM221A', 'MKRN3', 'MKRN1', 'MKRN9P', 'CCDC70', 'PSCA', 'DCUN1D2', 'LYZL6', 'CA8')
# gene_info <- bitr(gene_id, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
# which(!gene_id %in% gene_info$SYMBOL)
# gene_id[147]<-'CERT1'
# gene_id[212]<-'NIBAN1'
# #tox3
# gene_id<-c('KANK2', 'PRKAA2', 'GTSE1', 'FAM74A1', 'FAM74A4', 'FAM74A3', 'FAM74A7', 'FAM74A6', 'PGGHG', 'ADAM33', 'CDX1', 'SNAPC3', 'THPO', 'SPRED2', 'CBWD1', 'CBWD3', 'CBWD5', 'CBWD2', 'CBWD6', 'HEPACAM', 'UBXN2B', 'OIT3', 'G6PC2', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7D2', 'OR7G3', 'OR7G2', 'OR7G1', 'SFI1', 'PPP2CA', 'EIF2B2', 'UBA6', 'CAGE1', 'CEP19', 'NTN4', 'NEK6', 'ACD', 'KYNU', 'ZSCAN16', 'IST1', 'GLRX3', 'BYSL', 'CCNH', 'RB1', 'SPAG8', 'MCRS1', 'GATA2', 'GATA5', 'PNKP', 'TSTD2', 'EFHC1', 'BEND2', 'KDM1A', 'TLK2', 'SCEL', 'MRFAP1L1', 'MRFAP1', 'HSFY2', 'HSFY1', 'HSFY1P1', 'ESRRG', 'DRC3', 'NTM', 'OPCML', 'NDUFS1', 'NEIL2', 'HELLS', 'CCNC', 'CCDC103', 'BLOC1S5-TXNDC5', 'TXNDC5')
# #tox2
# gene_id<-c('C1orf109', 'CDKN2C', 'KLHL32', 'CIR1', 'GTSE1', 'POMP', 'UBTD2', 'DMTN', 'TAF1B', 'COG8', 'DTX2', 'CHCHD2', 'ZBED1', 'QARS', 'RB1', 'SEC23B', 'G2E3', 'SNAPC5', 'GSTT4', 'GSTTP2', 'CCDC102B', 'ZNF138', 'ZNF92', 'ZNF738', 'ZNF254', 'ZNF506', 'ZNF99', 'ZNF430', 'ZNF676', 'ZNF493', 'ZNF626', 'ZNF708', 'ZNF681', 'ZNF729', 'ZNF100', 'ZNF675', 'ZNF91', 'ZNF107', 'ZNF429', 'ZNF678', 'ZNF90', 'ZNF117', 'ERV3-1-ZNF117', 'ZNF730', 'ZNF695', 'ZNF718', 'ZNF66', 'ZNF732', 'KRTAP10-3', 'KRTAP10-9', 'KRTAP10-1', 'KRTAP10-11', 'KRTAP10-2', 'KRTAP10-4', 'KRTAP10-10', 'KRTAP10-12', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-5', 'CIB3', 'CPSF1', 'FAM74A1', 'FAM74A4', 'FAM74A3', 'FAM74A7', 'FAM74A6', 'ZNF572', 'PGGHG', 'LIMS3', 'LIMS4', 'LIMS3-LOC440895', 'LIMS1', 'GCM2', 'PDE4C', 'PDE4A', 'GNL3', 'TADA2B', 'NUFIP2', 'EMSY', 'XRN2', 'ACTN3', 'ACTN4', 'SPRED2', 'UBXN2B', 'OIT3', 'C22orf31', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G2', 'OR7G3', 'OR7D2', 'OR7E37P', 'OR7E5P', 'OR7E12P', 'OR7E14P', 'OR7E47P', 'OR7E2P', 'OR7G1', 'OR1M1', 'THPO', 'SLC25A30', 'XRCC4', 'SFI1', 'RCAN2', 'KIAA0825', 'SFN', 'RBCK1', 'SPIC', 'AGBL2', 'FCAR', 'FKBP6', 'PPP2CA', 'PAF1', 'RHOA', 'CAGE1', 'ZNF23', 'ZNF420', 'ZNF619', 'ZNF487', 'ZNF700', 'MYBL2', 'SIGLEC6', 'SIGLEC17P', 'CD33', 'ZNF202', 'TRIP10', 'RCOR3', 'BANP', 'SNX29P2', 'PPP2R1A', 'FAM89A', 'BLOC1S6', 'CKLF', 'CEP19', 'VPS28', 'TRIM74', 'TRIM73', 'TRIM50', 'GPANK1', 'CCNH', 'BYSL', 'GLRX3', 'KYNU', 'PDE4D', 'SNX8', 'HSPA8', 'MGC4859', 'NME7', 'GATA2', 'GATA5', 'GRN', 'HYDIN', 'HYDIN2', 'LILRB1', 'LILRA3', 'LILRA1', 'LILRA2', 'LILRB2', 'LILRB3', 'LILRB4', 'LILRB5', 'BEND2', 'MID2', 'HGS', 'PRMT7', 'PLCD4', 'XAB2', 'NLRC4', 'PACSIN2', 'PACSIN1', 'EML1', 'FAM129A', 'LCK', 'MYL6', 'MYL6B', 'RNH1', 'LDLRAP1', 'TLK2', 'MCMBP', 'CPNE6', 'ACD', 'ESRRG', 'DNAJB2', 'DRC3', 'GABPB1', 'NTM', 'OPCML', 'MORN3', 'ASPH', 'FXR2', 'HEXIM2', 'LOC105371795', 'NEIL2', 'CCNC', 'BLOC1S5-TXNDC5', 'TXNDC5', 'MKRN3', 'MKRN1', 'MKRN9P')
# gene_info <- bitr(gene_id, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
# which(!gene_id %in% gene_info$SYMBOL)
# gene_id[which(!gene_id %in% gene_info$SYMBOL)]
# gene_id[14]<-'QARS1'
# gene_id[177]<-'NIBAN1'
# #tox1
# gene_id<-c('C1orf109', 'CDKN2C', 'SHTN1', 'GTSE1', 'UBTD2', 'G2E3', 'SNAPC5', 'CCDC102B', 'KRTAP10-11', 'KRTAP10-9', 'KRTAP10-1', 'KRTAP10-2', 'KRTAP10-5', 'KRTAP10-4', 'KRTAP10-12', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-3', 'KRTAP10-10', 'KRTAP10-8', 'TADA2B', 'IGFN1', 'CLIP1', 'MON2', 'ZFP3', 'ZNF829', 'ZNF192P1', 'ZNF112', 'ZNF846', 'ZNF3', 'ZNF252P', 'ZNF582', 'ZNF723', 'ZFP37', 'ZNF568', 'ZNF891', 'ZNF594', 'ZNF84', 'ZNF735', 'ZNF727', 'ZNF527', 'ZNF883', 'ZNF506', 'ZNF354A', 'ZNF286A', 'ZNF286B', 'ZNF879', 'ZNF658', 'ZNF682', 'ZNF714', 'ZNF675', 'ZNF658B', 'ZNF540', 'ZNF624', 'ZNF436', 'ZNF660', 'TP53BP1', 'GRHL2', 'NUFIP2', 'MPP5', 'KLHL38', 'PLCB4', 'ACTN3', 'ACTN4', 'UBXN2B', 'OR7A17', 'OR7A10', 'OR7A5', 'OR7C1', 'OR7D4', 'OR7C2', 'OR7E24', 'OR7G2', 'OR7G3', 'OR7E37P', 'OR7E5P', 'OR7D2', 'CCZ1P-OR7E38P', 'OR7E12P', 'OR7E14P', 'OR7E47P', 'OR7E2P', 'OR7G1', 'OR1M1', 'HECTD3', 'PRKAG1', 'MALSU1', 'PPP2CA', 'HDAC7', 'PAF1', 'HSFY2', 'HSFY1', 'HSFY1P1', 'C1orf116', 'RHOA', 'RRS1', 'DEUP1', 'C1QTNF1', 'AGTR1', 'CNTD1', 'CPSF4', 'CD300LG', 'GFOD2', 'DHRS1', 'PRDX3', 'ZNF23', 'ZNF619', 'ZNF420', 'ZNF700', 'PHF11', 'SETDB2-PHF11', 'TEX2', 'ZFAND2B', 'MLX', 'TRIM74', 'TRIM73', 'TRIM50', 'MAPRE3', 'NT5C3B', 'CMBL', 'STARD7', 'TDP1', 'PBX3', 'TSNAX', 'TSNAX-DISC1', 'LILRB1', 'LILRA3', 'LILRA1', 'LILRA2', 'LILRB2', 'LILRB3', 'LILRB4', 'LILRB5', 'BEND2', 'RB1', 'DPP8', 'PLCD4', 'XAB2', 'ACOX3', 'RPN2', 'PROM1', 'SYT17', 'SLCO1C1', 'NFYB', 'TDRD7', 'PACSIN2', 'PACSIN1', 'PSMD2', 'TRIM32', 'EML1', 'FAM129A', 'LCK', 'CA8', 'RNH1', 'LDLRAP1', 'MCMBP', 'VCL', 'ACD', 'ESRRG', 'ARNTL', 'PARD6B', 'CATIP', 'NTM', 'OPCML', 'MORN3', 'ASPH', 'CNPY2', 'HEXIM2', 'LOC105371795', 'EMSY', 'NEIL2', 'CCNC', 'BLOC1S5-TXNDC5', 'TXNDC5', 'CRYBA1', 'MDM2')
# gene_info <- bitr(gene_id, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
# gene_id[which(!gene_id %in% gene_info$SYMBOL)]
# gene_id[151]<-'NIBAN1'



#tox8
gene_id<-c('TADA2B','MDM2','CKB','CKM','RHOA','UBTD2')
#tox7
gene_id<-c('EMSY','EIF2B2','CKB','CKM','ZMIZ2','GTSE1','ACD','TADA2B','AP4M1','NTM','OPCML','MORN3','HEPACAM','CCNC','MDM2')
#tox6
gene_id<-c('TADA2B','UBTD2','NUFIP2','RHOA','ZNF23','XRCC4','MLLT11','TXNDC5','BLOC1S5-TXNDC5')
#tox5
gene_id<-c('GTSE1','PRKAG1','MORN3','KANK2','HMBOX1','TPM3','TPM3P9','TPM4','DUS1L','ING3','HECTD3','MALSU1','PPP2CA','HDAC7','XRN2','WASHC3','SLA','MLX','PLCD4','EML1','ESRRG','EYA3','ZNF410','QARS','SNRNP25','PGGHG','SNAPC3','RRS1','NTM','OPCML','MAPRE3','STARD7','OTUB1','SNF8','GP1BA','DRC3','ARNTL','KIAA1586','PARD6B','SCAF1','RBP3','DTX2','COG4','GCM2','CLIP1','NUFIP2','UBXN2B','OR7A17','OR7A10','OR7A5','OR7C1','OR7D4','OR7C2','OR7E24','OR7G3','OR7D2','OR7G2','OR7E14P','OR7G1','C1orf116','SIGLEC6','SIGLEC17P','CD33','ZBTB16','ZNF202','TRIP10','RIMBP2','PRKACB','C1QTNF6','CMBL','PBX3','KYNU','INTS12','SPAG8','CPNE4','PNKP','NELFCD','LILRB1','LILRB2','LILRB3','LILRB4','LILRB5','LILRA1','LILRA2','LILRA3','VPS33A','CFAP70','MTMR4','PSMD2','ZSWIM2','PACSIN2','PACSIN1','TRIM32','TGM2','MCMBP','CPNE6','MKRN3','CRYBA1','CIB3','STRIP2','RNF181','MID2','CA8','PRMT7','DMTN','SEC23B','PCNP','KRTAP10-11','KRTAP10-9','KRTAP10-1','KRTAP10-2','KRTAP10-5','KRTAP10-4','KRTAP10-12','KRTAP10-6','KRTAP10-7','KRTAP10-3','KRTAP10-10','ACD','CCDC84','CPSF1','CCL16','BAGE2','BAGE3','BAGE5','BAGE4','BAGE','FAM74A1','FAM74A4','FAM74A3','FAM74A7','FAM74A6','RBBP7','NELL1','TADA2B','PHF12','UTP14C','UTP14A','ALG11','STAC2','CD22','DGKG','CDX1','SMCO1','MON2','ACTN3','ACTN4','SMARCB1','XRCC4','KIAA0825','FCAR','UEVLD','COPS8','CNPY2','TNFAIP8L2-SCNM1','SCNM1','PRKAG2','C1QTNF1','CD300LG','ZNF772','ZNF419 ','ZNF773','ZIK1','ZNF551','ZNF814','ZNF776','ZNF587','ZNF417','ZNF418','ZNF304','NAV1','KDM8','INPP5K','TUBG1','TUBG2','ARIH2','USB1','FAM89A','CEP19','RFPL3','RFPL2','RFPL1','RFPL3S','RFPL1S','XPO5','CIAPIN1','GLRX3','EME1','ASH2L','SLC22A17','PDE4D','MINDY1','GABRG1','NR1H2','SNPH','ZC3HC1','NELFA','SMAD1-AS1','SMAD1','GATA2','TSTD2','NOP58','CCNH','G2E3','RNH1','ASPH','NMRK1','FXR2','DCUN1D4','EMSY','CCDC103','FAM221A','CDH11','HGS')
#tox4
gene_id<-c('ACD', 'C1orf109', 'EML1', 'GYS1', 'LCK', 'RIMBP2', 'RNH1', 'SNAPC3', 'TPM3', 'TPM3P9', 'TPM4', 'VPS39', 'ALG11', 'ARFGAP1', 'ARPC3', 'CAGE1', 'CDCA7L', 'CDX1', 'CHCHD2', 'CIR1', 'COPB1', 'CPSF4', 'DCUN1D2', 'DGKG', 'EIF4EBP1', 'FAM221A', 'HMBOX1', 'KANK2', 'KLHL32', 'KRTAP10-8', 'KYAT1', 'KYNU', 'LYZL6', 'MBD4', 'MCM2', 'MKRN1', 'MKRN3', 'MKRN9P', 'NELL1', 'OR7A10', 'OR7A17', 'OR7A5', 'OR7C1', 'OR7C2', 'OR7D2', 'OR7D4', 'OR7E14P', 'OR7E24', 'OR7E2P', 'OR7E37P', 'OR7E47P', 'OR7E5P', 'OR7G1', 'OR7G2', 'OR7G3', 'PMS1', 'POPDC2', 'RASA1', 'RIPPLY1', 'RTN1', 'RWDD4', 'SHTN1', 'SMYD3', 'STAC', 'STAC2', 'THADA', 'UBA6', 'UBXN2B', 'WDR61', 'ZFYVE26', 'ZNF23', 'ZNF420', 'ZNF487', 'ZNF619', 'ZNF700', 'ZNF844', 'ZNF20', 'ZNF625-ZNF20', 'ZNF625', 'ZNF440', 'ZNF439', 'ZNF878', 'ZNF799', 'ZNF627', 'ZNF136', 'ZNF443', 'ZNF563', 'ZNF669', 'ZNF101', 'ZNF709', 'AHNAK', 'ATG12', 'ATG9A', 'BLOC1S5-TXNDC5', 'C11orf49', 'CA8', 'CARF', 'CCDC103', 'CCDC70', 'CCNH', 'CD33', 'CDC34', 'CINP', 'CKB', 'CKM', 'CNPY2', 'COL20A1', 'COL4A3BP', 'CPNE6', 'CPSF1', 'DMTN', 'FAM129A', 'G2E3', 'GPR68', 'HECTD3', 'ING3', 'KRTAP10-1', 'KRTAP10-10', 'KRTAP10-11', 'KRTAP10-12', 'KRTAP10-2', 'KRTAP10-3', 'KRTAP10-4', 'KRTAP10-5', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-9', 'LCN8', 'MRFAP1', 'MRFAP1L1', 'NMRK1', 'NUP62CL', 'OIT3', 'PDE4A', 'PDE4C', 'PGGHG', 'PNKP', 'PPP6R2', 'PSCA', 'RABL6', 'RBBP5', 'RRS1', 'SFI1', 'SIGLEC17P', 'SIGLEC6', 'SMAD1', 'SMAD1-AS1', 'SNPH', 'THAP10', 'TRAPPC2L', 'TXNDC5', 'UBTD1', 'USHBP1', 'UTP14C', 'UTP14A', 'YJU2', 'ZBTB16', 'ZNF202', 'ZNF572', 'C21orf59-TCP10L', 'ASH2L', 'BANP', 'C4orf33', 'CCDC198', 'CCNC', 'CFAP298', 'COQ3', 'CUTC', 'CYREN', 'DTX2', 'FAM74A1', 'FAM74A3', 'FAM74A4', 'FAM74A6', 'FAM74A7', 'GLRX3', 'GPANK1', 'HOXA1', 'KAZN', 'KIAA0825', 'LMBR1L', 'LYNX1-SLURP2', 'MCFD2', 'MCRS1', 'MID2', 'NAGK', 'NEK6', 'NTN4', 'PACSIN1', 'PACSIN2', 'PDE4D', 'POMP', 'PPP2CA', 'PRKAG1', 'RBX1', 'SLA', 'SNX29P2', 'TAF1B', 'TBC1D23', 'TSTD2', 'TTC14', 'ZBTB49','ZNF830')
#tox3
gene_id<-c('SNAPC3','OR7A17','OR7A10','OR7A5','OR7C1','OR7D4','OR7C2','OR7E24','OR7D2','OR7G3','OR7G2','OR7G1','THPO','G6PC2','PPP2CA','UBA6','RB1','SPAG8','DRC3','BLOC1S5-TXNDC5','TXNDC5','CDX1','SPRED2','SFI1','ACD','IST1','CCNH','PNKP','TSTD2','SCEL','PGGHG','CBWD1','CBWD3','CBWD5','CBWD2','CBWD6','HEPACAM','UBXN2B','NEK6','GLRX3','BYSL','BEND2','MRFAP1L1','MRFAP1','HSFY2','HSFY1','HSFY1P1','ESRRG','CCNC','CCDC103')
#tox2
gene_id<-c('CKLF','TADA2B','G2E3','FXR2','ASPH','MCMBP','PRMT7','VPS28','CEP19','FAM89A','ZNF23','ZNF420','ZNF619','ZNF487','ZNF700','RHOA','XRCC4','SPRED2','THPO','OR7A17','OR7A10','OR7A5','OR7C1','OR7D4','OR7C2','OR7E24','OR7G2','OR7G3','OR7D2','OR7E37P','OR7E5P','OR7E12P','OR7E14P','OR7E47P','OR7E2P','OR7G1','OR1M1','UBXN2B','NUFIP2','PDE4A','PDE4C','ZNF572','CIB3','CCDC102B','SNAPC5','RB1','UBTD2','CCNH','NTM','OPCML','DNAJB2','ACD','RNH1','MYL6','MYL6B','LCK','FAM129A','EML1','BEND2','LILRB1','LILRB2','LILRB3','LILRB4','LILRB5','LILRA3','LILRA1','LILRA2','NME7','KYNU','GLRX3','BLOC1S6','ZNF202','SIGLEC6','SIGLEC17P','CD33','CAGE1','PAF1','PPP2CA','RBCK1','SFN','KIAA0825','GNL3','LIMS3','LIMS4','LIMS3-LOC440895','LIMS1','CPSF1','KLHL32','RCAN2','MKRN3','MKRN1','MKRN9P','TXNDC5','BLOC1S5-TXNDC5','CCNC','NEIL2','MORN3','GABPB1','DRC3','ESRRG','CPNE6','PACSIN1','PACSIN2','NLRC4','XAB2','MID2','HYDIN','HYDIN2','GRN','BYSL','PPP2R1A','SNX29P2','BANP','RCOR3','TRIP10','MYBL2','FKBP6','FCAR','SPIC','SFI1','SLC25A30','C22orf31','OIT3','ACTN3','ACTN4','XRN2','EMSY','GCM2','KRTAP10-1','KRTAP10-2','KRTAP10-3','KRTAP10-4','KRTAP10-5','KRTAP10-6','KRTAP10-7','KRTAP10-9','KRTAP10-10','KRTAP10-11','KRTAP10-12','ZNF138','ZNF92','ZNF738','ZNF254','ZNF506','ZNF99','ZNF430','ZNF676','ZNF493','ZNF626','ZNF708','ZNF681','ZNF729','ZNF100','ZNF675','ZNF91','ZNF107','ZNF429','ZNF678','ZNF90','ZNF117','ERV3-1-ZNF117','ZNF730','ZNF695','ZNF718','ZNF66','ZNF732','GSTT4','GSTTP2','ZBED1','CHCHD2','TAF1B','POMP','CIR1','CDKN2C','C1orf109','SNX8','GPANK1','TRIM74','TRIM73','TRIM50','FAM74A1','FAM74A4','FAM74A3','FAM74A7','FAM74A6','SEC23B','QARS','DTX2','COG8','DMTN')
#tox1
gene_id<-c('EMSY','MCMBP','TADA2B','UBTD2','ACD','RHOA','CLIP1','G2E3','GTSE1','EML1','NFYB','HDAC7','PPP2CA','SNAPC5','HSFY1','HSFY2','HSFY1P1','CCNC','CNPY2','MORN3','ARNTL','ESRRG','VCL','LDLRAP1','RNH1','CA8','PACSIN2','PACSIN1','DPP8','TDP1','ZFAND2B','PHF11','SETDB2-PHF11','ZNF23','ZNF619','ZNF420','ZNF700','RRS1','MALSU1','NUFIP2','ZFP3','ZNF829','ZNF192P1','ZNF112','ZNF846','ZNF3','ZNF252P','ZNF582','ZNF723','ZFP37','ZNF568','ZNF891','ZNF594','ZNF84','ZNF735','ZNF727','ZNF527','ZNF883','ZNF506','ZNF354A','ZNF286A','ZNF286B','ZNF879','ZNF658','ZNF682','ZNF714','ZNF675','ZNF658B','ZNF540','ZNF624','ZNF436','ZNF660','KRTAP10-1','KRTAP10-2','KRTAP10-3','KRTAP10-4','KRTAP10-5','KRTAP10-6','KRTAP10-7','KRTAP10-8','KRTAP10-9','KRTAP10-10','KRTAP10-11','KRTAP10-12','NEIL2','CATIP','OR7A17','OR7A10','OR7A5','OR7C1','OR7D4','OR7C2','OR7E24','OR7G2','OR7G3','OR7E37P','OR7E5P','OR7D2','CCZ1P-OR7E38P','OR7E12P','OR7E14P','OR7E47P','OR7E2P','OR7G1','OR1M1','PLCD4','RB1','BEND2','TSNAX','TSNAX-DISC1','PBX3','STARD7','MLX','TEX2','C1QTNF1','C1orf116','HECTD3','ACTN3','ACTN4','TP53BP1','IGFN1','CCDC102B','SHTN1','CDKN2C','C1orf109')


gene_info <- bitr(gene_id, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
which(!gene_id %in% gene_info$SYMBOL)
gene_id[which(!gene_id %in% gene_info$SYMBOL)]
target_gene_id <- as.character(gene_info[, 2])
length(target_gene_id)

###画图函数，这里面是GO富集分析的一些常见图形
plot_fun <- function(ego, prefix){
  png(paste(prefix, "barplot.png", sep = "_"))
  barplot(ego, showCategory = 30)
  dev.off()
  png(paste(prefix, "dotplot.png", sep = "_"))
  dotplot(ego, showCategory = 30)
  dev.off()
  png(paste(prefix, "plotGOgraph.png", sep = "_"))
  plotGOgraph(ego)
  dev.off()
  png(paste(prefix, "goplot.png", sep = "_"))
  goplot(ego)
  dev.off()
  png(paste(prefix, "emapplot.png", sep = "_"))
  emapplot(ego, showCategory = 30)
  dev.off()
}

###GO富集分析
ego_MF <- enrichGO(OrgDb = "org.Hs.eg.db",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "MF",
                   readable = TRUE)
ego_result_MF <- as.data.frame(ego_MF)

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "CC",
                   readable = TRUE)
ego_result_CC <- as.data.frame(ego_CC )

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "BP",
                   readable = TRUE)
ego_result_BP <- as.data.frame(ego_BP)

# write.table(ego_result_MF, "tox1-go_enrich_MF.xls", quote = F, sep = "\t", col.names = T, row.names = F)
# ego_result_MF <- ego_result_MF[1:10,]
# write.table(ego_result_BP, "FL-go_enrich_BP.xls", quote = F, sep = "\t", col.names = T, row.names = F)
# ego_result_BP <- ego_result_BP[1:10,]
# write.table(ego_result_CC, "tox3-go_enrich_CC.xls", quote = F, sep = "\t", col.names = T, row.names = F)
# ego_result_CC <- ego_result_CC[1:10,]


go_enrich_df <- data.frame(ID=c(as.character(ego_result_BP$ID), as.character(ego_result_CC$ID), as.character(ego_result_MF$ID)),
                           Description=c(as.character(ego_result_BP$Description), as.character(ego_result_CC$Description), as.character(ego_result_MF$Description)),
                           Qval = c(ego_result_BP$qvalue, ego_result_CC$qvalue, ego_result_MF$qvalue),
                           GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", nrow(ego_result_BP)), rep("cellular component", nrow(ego_result_CC)),
                                         rep("molecular function", nrow(ego_result_MF))), levels=c("molecular function", "cellular component", "biological process")))


###绘制GO富集柱状图
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=9, n_char=60){  
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char)){    
    if (nchar(x) > n_char){
      x <- substr(x, 1, n_char)
      x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                 collapse=" "), "...", sep="")    
      return(x)
    }
  }
  else{
      return(x)
    }
}

labels=(sapply(go_enrich_df$Description,
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))
#导入字体
#extrafont::loadfonts()
library(ggplot2)
go_enrich_df$padj <- ifelse(go_enrich_df$Qval < 0.05, '*', '')
p <- ggplot(data = go_enrich_df, aes(x = number, y = -log(Qval), fill = type)) +
  scale_fill_manual(values=c("cellular component"= "#F8766D", "biological process"= "#00BFC4", "molecular function"= "#00BA38"))+
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  geom_text(aes(label = padj), vjust=0.75, hjust = -0.25, color = 1) +
  #scale_fill_manual(values = CPCOLS) + 
  theme_bw() + 
  scale_x_discrete(labels=labels) +
#  theme(axis.text=element_text(face = "bold", color="gray50")) +
  theme(text=element_text(family="Arial"))+
  labs(title = "Ece-VI")+
  xlab(" ")+
  scale_y_continuous(name = "-log10(q-value)",breaks=seq(0,30,4))+
  theme(panel.grid.major=element_line(colour=NA),           
        panel.grid.minor=element_blank())+                  #移除网格线
  theme(panel.background = element_rect(colour = "black",size = 1.5))+    #边框加粗
  theme(plot.title = element_text(size = 16),
        text = element_text(size = 15),
       axis.title = element_text(size = 15),
        axis.text.x=element_text(size = 15))
p
ggsave(filename = "F1A-tox6.pdf", plot = p, height=3.62, width=10,dpi = 300)#5.5/9..2/8..
ggsave(filename = "tox-go_enrich_of_diffgene.png", type="cairo-png", plot=p, height = 6, width = 8)

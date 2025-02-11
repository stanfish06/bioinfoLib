import numpy as np
import scanpy as sc

# Idse's gene list: 3/22/24
coremarkers = ['SOX2','TBXT','SOX17','ISL1']
epiblast = ['PODXL','SOX2','NANOG','POU5F1','KLF4','OTX2','ESRG','DPPA4','USP44']
streak = ['TBXT','MIXL1','EOMES']
mesoderm = ['GSC','TBX6','MESP1','GATA6','GATA4','MSGN1','PDGFRA','EPHA4','ZIC3','LHX1','DLL3','RSPO3','FOXC1'] # ,'EPHA4','ZIC3' Zhou..Lanner..Chien embryonic mesoderm
endoderm = ['OTX2','PRDM1','SOX17','FOXA2','HHEX','LHX1','TTR','GATA4','KIT','CXCR4']
ectoderm = ['SOX1', 'SOX2', 'SOX3','NES']
TE = ['CGA','XAGE3','PGF','WNT3A']
amnion = ['ISL1','BMP4','GABRP','HEY1','LYPD1','WNT6']#,'WFDC2','ANXA3','CDO1'] #Zhou..Lanner...
sexmeso = ['PITX1','TBX4','CYB5A','MLLT3','ALCAM'] #Zhou..Lanner...
ysmeso = ['ETS1','ANXA1','POSTN','DCN']
exMC = ['VIM','COL3A1','HGF','COL6A2','HAND2','COL1A1','POSTN']
extraembryonic = ['CDX2', 'GATA2','GATA3','HAND1','TBX3','TFAP2A','DLX5','KRT18','KRT7','TP63'] + TE + amnion + sexmeso + ysmeso + exMC;
PGC = ['PRDM1','SOX17','TFAP2C','NANOG','ALPP','DPPA4','NANOS3','LAMA4','KLF4','KIT','CXCR4','DDX4','DAZL','DPPA3'] # DDX4=VASA, DPPA3=stella
endothelium = ['PECAM1', 'MEF2C']
blood = ['RUNX1','GATA1','HBE1']
heart = ['NKX2-5','MESP1','TBX1','TBX5','ISL1','TNNT2']
lung = ['NKX2-1']
organs = endothelium + blood + heart + lung
HOX = ['HOX[ABCD][1-9]+'];
germmarkers =  list(np.unique(epiblast + mesoderm + endoderm + ectoderm + extraembryonic + PGC))
allmarkers = list(np.unique(germmarkers + organs + HOX))
coreligands = ['DKK[1-9]','LEFTY[1-2]','NODAL','WNT3','WNT3A','WNT6','WNT5A','WNT5B','WNT8A','BMP2','BMP4','BMP7','CER1','NOG','FGF2','FGF4','FGF8','FGF17']
TGFbtransduction = ['SMAD[1-9]','ACVR[1-2][ABC]*','BMPR[1-2][ABC]*','TGFBR','TDGF1']
wnt = ['WNT[1-9]+[AB]*','RSPO[1-9]','FZD[1-9]*','DVL']
othersignaling = ['EGF','EGFR']
FGFR = ['FGFR[1-4]']
extracell_FGFact = ['KL','GLG1','FLRT', 'KAL1'] # klotho, cfr, flrt, anosmin-1
extracell_FGFinh = ['FGFRL1','IL17RD'] # FGFRL1, SEF
intracell_ERKinh = ['SPRY[1-4]','DUSP[1-9]','DUSP10','DUSP14','DUSP26','PEBP1']
intracell_ERKact = ['MAP2K[1-2]','MAPK3','MAPK1','ARAF','BRAF','RAF1'] #K3 = ERK1, #K1 = ERK2, MEK1/2 are phosporylated by RAS
HSPG_synthesis = ['EXT1','EXT2','UGDH'] # EXT catalyze HS side chain elongation, UGDH : Garcia-Garcia & Anderson, makes side chains
HSPG = ['HSPG2', 'AGRIN', 'SDC[1-4]', 'GPC[1-9]'] # HSPG2: perlecan (secreted), agrin=secreted, syndecan=transmembrane, glypican (GPI-anchored)
HSPG_desulf = ['SULF'] # 6-O sulfatases - inhibit FGF
HSPG_shedding = ['HPSE', 'HTRA1'] # heparanase
HSPG_sulfotransferase = ['HS6ST','HS2ST', 'HS3ST', 'NDST']
FGF_inhibition = HSPG_desulf + extracell_FGFinh + intracell_ERKinh
FGF_activation = FGFR + extracell_FGFact + HSPG_synthesis + HSPG + HSPG_shedding + HSPG_sulfotransferase + intracell_ERKact
FGF_modulation = FGF_inhibition + FGF_activation
allsignaling = coreligands + TGFbtransduction + FGF_modulation + wnt + othersignaling;
EMT = ['CDH[1-2]','VIM','EPCAM','SNAI[1-2]','PAI[1-3]']
allgenes = list(np.unique(allmarkers + EMT + allsignaling))
allgenes_re_symbol = list(np.unique(allmarkers + EMT + allsignaling))
allgenes = ['^' + g + '($| )' for g in allgenes]

def normalize_and_select_hvg(
    adata,
    target_sum=1e4,
    n_top_genes=2000,
):
    adata.layers["count"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw=  adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,subset=True,flavor="seurat_v3",layer="count")

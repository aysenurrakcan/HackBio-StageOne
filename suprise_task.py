# Figure a

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = """gene,HBR_1,HBR_2,HBR_3,UHR_1,UHR_2,UHR_3
SULT4A1,375,343.6,339.4,3.5,6.9,2.6
MPPED1,157.8,158.4,162.6,0.7,3,2.6
PRAME,0,0,0,568.9,467.3,519.2
IGLC2,0,0,0,488.6,498,457.5
IGLC3,0,0,0,809.7,313.8,688
CDC45,2.6,1,0,155,152.5,149.9
CLDN5,77.6,88.5,67.2,1.4,2,0
PCAT14,0,0,1.2,139.8,154.4,155.1
RP5-1119A7.17,53,57.6,51.9,0,0,0
MYO18B,0,0,0,59.5,84.2,56.5
RP3-323A16.1,0,0,1.2,51.9,76.2,53.1
CACNG2,42.7,35,56.6,0,1,0
"""

from io import StringIO
df = pd.read_csv(StringIO(data), index_col=0)

sns.set(font_scale=0.9)
g = sns.clustermap(
    df,
    cmap="Blues",      
    linewidths=0.5,       
    figsize=(6,5),
    cbar_kws={"label": "Normalized Expression"},
)

g.fig.suptitle("a. Clustered Heatmap of Differentially Expressed Genes", y=1.05)
plt.show()



 #Figure b
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_volcano_plot(data_url: str, lfc_thresh: float):
    """
    Generates a Volcano Plot from DEG data at the specified URL.

    data_url: The URL of the dataset.
    lfc_thresh: The threshold value for log2FoldChange.
    """

    C_MAP = {'up': 'green', 'down': 'orange', 'ns': 'grey'}
    LFC_C, P_ADJ_C, SIG_C, LOG_P_C = 'log2FoldChange', 'PAdj', 'significance', '-log10(Padj)'

    try:
        df = pd.read_csv(data_url)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

  
    df = (df
          .assign(
              **{P_ADJ_C: lambda x: x[P_ADJ_C].replace(0, x[P_ADJ_C][x[P_ADJ_C] > 0].min() or 1e-300)}
          )
          .assign(
              **{LOG_P_C: lambda x: -np.log10(x[P_ADJ_C])}
          )
          .sort_values(by=SIG_C, ascending=False)
    )

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 8))

    for cat, color in C_MAP.items():
        sub = df[df[SIG_C] == cat]

        ax.scatter(
            sub[LFC_C],
            sub[LOG_P_C],
            s=30,
            alpha=0.8,
            color=color,
            label={'up': 'Upregulated', 'down': 'Downregulated', 'ns': 'Not Significant'}[cat],
            edgecolors='none'
        )

    ax.axvline(x=lfc_thresh, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(x=-lfc_thresh, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

    ax.set_xlabel('log2FoldChange', fontsize=14, labelpad=10)
    ax.set_ylabel(r'$-\log_{10}(Padj)$', fontsize=14, labelpad=10)
    ax.set_title('Volcano Plot of Differential Gene Expression', fontsize=16)

    ax.set_xlim(df[LFC_C].min() - 1, df[LFC_C].max() + 1)
    ax.set_ylim(0, df[LOG_P_C].max() + 10)

    ax.legend(title='Significance', loc='upper right', frameon=True, fontsize=10)

    plt.show()

DATA_URL = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
LFC_THRESHOLD = 1
generate_volcano_plot(DATA_URL, LFC_THRESHOLD)


#Figure c
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_scatter_plot(data_url: str):
  
   
    RADIUS_COL, TEXTURE_COL, DIAGNOSIS_COL = 'radius_mean', 'texture_mean', 'diagnosis'
    
    COLOR_MAP = {
        'M': 'C0', 
        'B': 'C1'   
    }
    
    LABEL_MAP = {
        'M': 'Malignant (M)',
        'B': 'Benign (B)'
    }

    try:
        df = pd.read_csv(data_url)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    if not all(col in df.columns for col in [RADIUS_COL, TEXTURE_COL, DIAGNOSIS_COL]):
        print("Required columns (radius_mean, texture_mean, diagnosis) not found in the dataset.")
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 6))

    for diag_type, color in COLOR_MAP.items():
        subset = df[df[DIAGNOSIS_COL] == diag_type]
        
        ax.scatter(
            subset[RADIUS_COL],
            subset[TEXTURE_COL],
            s=40,
            alpha=0.8,
            color=color,
            label=LABEL_MAP[diag_type],
            edgecolors='none'
        )

    ax.set_xlabel(RADIUS_COL, fontsize=14, labelpad=10)
    ax.set_ylabel(TEXTURE_COL, fontsize=14, labelpad=10)
    ax.set_title('Texture Mean vs. Radius Mean by Diagnosis', fontsize=16)

    ax.legend(title='Diagnosis', loc='upper right', frameon=True, fontsize=10)
    
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)
    
    plt.show()

DATA_URL = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
generate_scatter_plot(DATA_URL)


#Figure D
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def generate_correlation_heatmap(data_url: str):
   
    FEATURES = [
        'radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 
        'smoothness_mean', 'compactness_mean'
    ]

    try:
        df = pd.read_csv(data_url)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    df_subset = df[FEATURES]
    
    correlation_matrix = df_subset.corr()

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(
        correlation_matrix,
        annot=True,             
        fmt=".1f",               
        cmap="Blues",            
        linewidths=.5,           
        linecolor='white',
        cbar_kws={'label': 'Correlation Coefficient'}, 
        ax=ax
    )

    ax.set_title('Correlation Heatmap of Key Mean Features', fontsize=16, pad=20)
    
    ax.tick_params(axis='x', rotation=90)
    ax.tick_params(axis='y', rotation=0)
    
    plt.tight_layout()
    plt.show()

DATA_URL = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
generate_correlation_heatmap(DATA_URL)



#Figure E
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

def generate_smoothness_compactness_scatter(data_url: str):
   
    SMOOTHNESS_COL = 'smoothness_mean'
    COMPACTNESS_COL = 'compactness_mean'
    DIAGNOSIS_COL = 'diagnosis'
    
    COLOR_MAP = {
        'M': 'C0',  
        'B': 'C1'   
    }
    
    LABEL_MAP = {
        'M': 'Malignant (M)',
        'B': 'Benign (B)'
    }

    try:
        df = pd.read_csv(data_url)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    required_cols = [SMOOTHNESS_COL, COMPACTNESS_COL, DIAGNOSIS_COL]
    if not all(col in df.columns for col in required_cols):
        print(f"Required columns {required_cols} not found in the dataset.")
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 6))

    for diag_type, color in COLOR_MAP.items():
        subset = df[df[DIAGNOSIS_COL] == diag_type]
        
        ax.scatter(
            subset[SMOOTHNESS_COL],
            subset[COMPACTNESS_COL],
            s=40,
            alpha=0.8,
            color=color,
            label=LABEL_MAP[diag_type],
            edgecolors='none'
        )

    ax.set_xlabel(SMOOTHNESS_COL, fontsize=14, labelpad=10)
    ax.set_ylabel(COMPACTNESS_COL, fontsize=14, labelpad=10)
    ax.set_title('Compactness Mean vs. Smoothness Mean by Diagnosis', fontsize=16)

    ax.legend(title='Diagnosis', loc='upper left', frameon=True, fontsize=10)
    
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)
    
    plt.show()

DATA_URL = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
generate_smoothness_compactness_scatter(DATA_URL)


#Figure F

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def generate_density_plot(data_url: str):
  
    
    AREA_COL = 'area_mean'
    DIAGNOSIS_COL = 'diagnosis'
    
    LABEL_MAP = {
        'M': 'Malignant (M)',
        'B': 'Benign (B)'
    }

    try:
        df = pd.read_csv(data_url)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    required_cols = [AREA_COL, DIAGNOSIS_COL]
    if not all(col in df.columns for col in required_cols):
        print(f"Required columns {required_cols} not found in the dataset.")
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sns.kdeplot(
        data=df,
        x=AREA_COL,
        hue=DIAGNOSIS_COL,
        hue_order=['M', 'B'], 
        fill=True,          
        alpha=0.2,           
        linewidth=2,
        common_norm=False,   
        palette=['C0', 'C1'], 
        ax=ax
    )
    
    legend = ax.get_legend()
    if legend:
        for text in legend.get_texts():
            current_label = text.get_text()
            if current_label in LABEL_MAP:
                text.set_text(LABEL_MAP[current_label])

    ax.set_xlabel(AREA_COL, fontsize=14, labelpad=10)
    ax.set_ylabel('Density', fontsize=14, labelpad=10)
    ax.set_title('Kernel Density Estimate of Area Mean by Diagnosis', fontsize=16)

    ax.legend_.set_title('Diagnosis')
    
    plt.tight_layout()
    plt.show()

DATA_URL = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
generate_density_plot(DATA_URL)

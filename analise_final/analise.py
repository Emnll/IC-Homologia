#%%
import pandas as pd
from collections import Counter

#%%
df = pd.read_csv("message (3).txt", sep='\t')
# KOG, KOG_OU_FUNC,         
#%%
if 'GO_ID' not in df.columns:
   print("A coluna 'GO_ID' não foi encontrada no arquivo.")


            # Handle rows where 'GO_ID' might be missing or empty
go_ids = df['GO_ID'].dropna()

            # Split the comma-separated values and flatten the list
all_go_ids = []
for ids in go_ids:
    all_go_ids.extend([go.strip() for go in ids.split(',')])

            # Count the frequency of each GO_ID
go_id_counts = Counter(all_go_ids)
print(go_id_counts )

            # Find the most common GO_ID
if go_id_counts:
    most_common_id, count = go_id_counts.most_common(1)[0]
    print(f"O GO_ID mais frequente é '{most_common_id}', com {count} ocorrências.")
# %%

import muspan as ms
import matplotlib.pyplot as plt
import numpy as np


# load in the dataset
example_domain=ms.datasets.load_example_domain('Xenium-Healthy-Colon')

# print the labels 
example_domain.print_labels()


# %%

example_domain.estimate_boundary(method='convex hull')

# visualise the dataset to make sure it's the right one
ms.visualise.visualise(example_domain,color_by='Cluster ID',marker_size=10,objects_to_plot=('collection','Cell boundaries'),show_boundary=True)


# %%
example_domain.print_labels('Cell ID')


# %% Querying everything associated with a cell ID

plt.close('all')

cell_id_of_interest = 'dfbfdjho-1'

short_hand_query_cell_id=('Cell ID',cell_id_of_interest) # clever stuff in backend for short hand

full_cell_id_query=ms.query.query(example_domain,('label','Cell ID'),'is',cell_id_of_interest)  # short hand will just make this query

ms.visualise.visualise(example_domain,marker_size=10,color_by='Cluster ID',objects_to_plot=full_cell_id_query)



# %% Querying everything associated with a cell ID and SOX9

plt.close('all')

query_cell_id_and_sox9=ms.query.query_container(short_hand_query_cell_id,'AND',('Transcript ID', 'Sox9'), domain=example_domain)

ms.visualise.visualise(example_domain,marker_size=10,color_by='Transcript ID',objects_to_plot=query_cell_id_and_sox9)


# %% Get all cells within 20µm distance of the cell with ID 'dfbfdjho-1'. Add these cells to a collection call 'Close to dfbfdjho-1'

plt.close('all')

# all cells
query_all_cells = ('collection','Cell boundaries')

# distance query - all cells (object 1) from cells with ID (object 2)
query_all_close_to_cell = ms.query.query(example_domain,('distance', ('boundary', query_all_cells, short_hand_query_cell_id)), '<' ,20)

ms.visualise.visualise(example_domain,marker_size=10,objects_to_plot=query_all_close_to_cell)


# %% add these to a collection 'Close to dfbfdjho-1'

example_domain.add_objects_to_collection(add_collection_to=query_all_close_to_cell,collection_name='Close to dfbfdjho-1')

print(example_domain)

ms.visualise.visualise(example_domain,marker_size=10,color_by='Cell ID',objects_to_plot=('collection','Close to dfbfdjho-1'))



# %% - EXTENSION PROBLEM : Get all Transcripts within the collection 'Close to dfbfdjho-1'.

in_collection_query=ms.query.query(example_domain,('collection',),'is','Close to dfbfdjho-1')

# first we need all the cell IDs in associated with the cells in this collection
all_ids,ids_objects = ms.query.get_labels(domain=example_domain, label_name='Cell ID')

# actually getting the object IDs of all objects in the collection
in_coll_object_ids=ms.query.interpret_query(in_collection_query)

_,_,indices_in_labels=np.intersect1d(in_coll_object_ids,ids_objects,return_indices=True)

these_cell_ids = all_ids[indices_in_labels]

print(these_cell_ids)

# %%  now we have the cell ids within this collection we can query all objects that are associated with these IDs

query_in_list_of_cell_ids=ms.query.query(example_domain,('label','Cell ID'),'in',these_cell_ids)

query_transcripts_associated_collection=ms.query.query_container(query_in_list_of_cell_ids,'AND',('collection', 'Transcripts'), domain=example_domain)

ms.visualise.visualise(example_domain,color_by='Transcript ID',marker_size=10,objects_to_plot=query_transcripts_associated_collection)



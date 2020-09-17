from collections import defaultdict

import numpy as np
import pandas as pd
from pastml.tree import read_tree, remove_certain_leaves
import random

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()


    parser.add_argument('--input_tree', default='C:/Users/Andrew/PycharmProjects/masters_project/newyork_subsample/collapsed_rooted_26AUG2020_DATED.nwk', type=str) # you need to put the path to your own tree for input and output
    parser.add_argument('--input_locs', default='C:/Users/Andrew/PycharmProjects/masters_project/tree_meta_org/group_metadata_0.tab', type=str) # metadata which includes IDs and locations
    parser.add_argument('--input_counts', default='C:/Users/Andrew/PycharmProjects/masters_project/25April2020_pergroup.csv', type=str)  # state case counts since April 25, 2020
    parser.add_argument('--size', default=200, type=int, help='number of tips in the subsampled tree (should be inferior to the one in the initial tree)') # this number must be smaller than the total number of sequences
    parser.add_argument('--output_ids', type=str, nargs='+', default=['C:/Users/Andrew/PycharmProjects/masters_project/newyork_subsample/ids1_DATED_SUB_200.nwk', 'C:/Users/Andrew/PycharmProjects/masters_project/tree0_subsample/ids2.nwk'])
    parser.add_argument('--output_trees', type=str, nargs='+', default=['C:/Users/Andrew/PycharmProjects/masters_project/newyork_subsample/tree1_DATED_SUB_200.nwk', 'C:/Users/Andrew/PycharmProjects/masters_project/tree0_subsample/tree2.nwk'])
    parser.add_argument('--output_stats', required=False, type=str, default='C:/Users/Andrew/PycharmProjects/masters_project/newyork_subsample/stats_DATED_SUB_200.table')
    params = parser.parse_args()

    # phylogenetic tree read from the newick file
    # Check out ETE3 tutorial for tree manipulation:
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html
    tree = read_tree(params.input_tree)


    # check out the tree.traverse() function, test printing node names etc.
    # preorder visits the parent nodes before their children, postorder visits children before their parent
    # for n in tree.traverse('preorder'):
    #   print(n.name, n.dist)

    # reading the metadata table
    loc_df = pd.read_csv(params.input_locs, index_col=0, sep='\t')
    loc_df.index = loc_df.index.map(str)

    # remove tree leaves with unknown states
    tree = remove_certain_leaves(tree, to_remove=lambda _: (_.name not in loc_df.index) or pd.isna(loc_df.loc[_.name, 'county']))

    ids = {_.name for _ in tree}
    # keeping only the rows of the table that correspond to tree leaves
    loc_df = loc_df.loc[loc_df.index.isin(ids), :]

    # group the table rows by state and count the number of sampled sequences per state
    # counting number of sequences in each state

    sampled_case_df = loc_df[['county', 'country']].groupby(['county']).count()
    sampled_case_df.columns = ['sampled_cases']
    print(sampled_case_df.head()) #ONLY ONE CASE IN NJ!! Despite the second biggest case number
    print(sampled_case_df.sum(axis=0, skipna=True))

    # here we are going to keep equal proportions of different states,
    # you would need to keep the proportions of states that you define according to declared cases
    # or some other external knowledge that you find pertinent

    # dataframes in python are similar to the R ones:
    # check out https://pandas.pydata.org/docs/getting_started/comparison/comparison_with_r.html

    #state_counts = open(params.input_counts)
    #state_counts = reader(state_counts)
    #state_counts = list(state_counts)

    #READING THE STATE CASE DATA
    case_df = pd.read_csv(params.input_counts, index_col = 0)
    sum_cases = case_df.sum(axis=0, skipna=True)

    #case_df = pd.DataFrame(index=sampled_case_df.index)
    left = params.size
    n = len(case_df) #I'm not sure if this is right
    avg_to_take = max(left / 50, 1) # I do not understand this, does this have to be num of states?
    print('Aiming to take {} samples per state.'.format(avg_to_take))
    # this is the target proportions for each state (here equal)

    #calculate the frequency by number of cases in state divided by total number of cases
    case_df['frequencies'] = case_df['cases'] / 1218666
    # the number of sequences we want to keep for each state - size multiplied by the frequency
    case_df['rescaled_cases'] = case_df['frequencies'] * params.size #THIS IS WHAT MAKE SENSE TO ME

    print(case_df[0:5])

    # the number of available sequences per each state from tree data. Merged.
    case_df = pd.merge(case_df, sampled_case_df, on = 'county') #NOT SURE IF THIS IS AN OKAY TO DO THIS


    #zcase_df['sampled_cases'] = sampled_case_df.loc[case_df.index, 'April25_cases']


    case_df.sort_values(by=['sampled_cases'], inplace=True, ascending=True)

    #print(case_df[0:5])


    # if there are not enough sequences sampled (less than we want to keep),
    # we are going to take a bit more of other state sequences to compensate for it
    remaining = {} # new dictionary to add states that have sequences remaining


    case_df['subsampled_cases'] = 0

    def subsample(size, freq_df):

        temp_df = freq_df[['frequencies', 'sampled_cases']]
        temp_df['rescaled_cases'] = temp_df['frequencies'] * size

        for i, (county, row) in enumerate(temp_df.iterrows()):
            avail = row['sampled_cases']
            #print(round(row['rescaled_cases']))
            like_to_take = np.round(row['rescaled_cases'])
            if like_to_take < 3:
                like_to_take = 3
            can_take = min(like_to_take,avail)
            case_df.loc[county,'subsampled_cases'] += can_take
            size -= can_take
            print(can_take, size, row['sampled_cases'])
            temp_df.loc[county, 'sampled_cases'] -= can_take
            #row['sampled_cases'] = row['sampled_cases'] - can_take
            #print(row['sampled_cases'])
        #print(temp_df['sampled_cases']) #this doesn't reflect the same as the row above
        not_empty = temp_df['sampled_cases'] > 0
        #print(not_empty)
        temp_df = temp_df[not_empty]
        #print(temp_df)
        if size > 0 and len(temp_df) > 0:
            subsample(size, temp_df)


    z = subsample(params.size, case_df)

    case_df.sort_values(by=['subsampled_cases'], inplace=True, ascending=True)
    # save the subsampling stats table
    case_df.to_csv(params.output_stats, sep='\t', index_label='index')

    total_subsampled = case_df.sum(axis=0, skipna=True)
    print(total_subsampled)

    print('There are a total of {} sampled counties.'.format(len(case_df)))
    print(case_df)
    print(case_df['subsampled_cases'].sum())

    print(case_df.columns)
    #case_df.loc["abroad",:] = ["AB", 0,0,0,0]


    # now let's select the ids of the sequences to keep and produce subsampled trees
    for (output_tree, output_ids) in zip(params.output_trees, params.output_ids):
        tree = read_tree(params.input_tree)
        # remove tree leaves with unknown states
        tree = remove_certain_leaves(tree, to_remove=lambda _: (_.name not in loc_df.index) or pd.isna(loc_df.loc[_.name, 'county']))

        # let's calculate how many sequences we need to remove for each country
        # and keep updating this number while pruning the tree
        case_df['remove_cases'] = case_df['sampled_cases'] - case_df['subsampled_cases']

        print(case_df)

        bin2tips = defaultdict(list)
        # go through the tree leaves, and group them by their branch length
        for _ in tree:
            # if there we can still remove sequences of this state, put this leaf in its branch length bin
            if case_df.loc[loc_df.loc[_.name, 'county'], 'remove_cases']:
                bin2tips[_.dist].append(_)
        # while there are still sequences to be removed let's keep pruning the tree
        while case_df['remove_cases'].sum():
            # pick a leaf with the minimal branch length (randomly if there are several such leafs)
            min_bin = min(bin2tips.keys())
            min_tips = bin2tips[min_bin]
            tip = np.random.choice(min_tips, 1)[0]
            min_tips.remove(tip)
            # if nothing is left for this bin, remove it
            if not min_tips:
                del bin2tips[min_bin]
            # if we can remove this sequence let's do it
            if case_df.loc[loc_df.loc[tip.name, 'county'], 'remove_cases'] > 0:
                case_df.loc[loc_df.loc[tip.name, 'county'], 'remove_cases'] -= 1
                # prune the leaf from the tree
                parent = tip.up
                parent.remove_child(tip)
                # if the leaf we just pruned had only one sibling,
                # we need to merge it with its parent and update its branch length accordingly
                if len(parent.children) == 1 and not parent.is_root():
                    grandparent = parent.up
                    child = parent.remove_child(parent.children[0])
                    if child.is_leaf() and child.dist in bin2tips and child in bin2tips[child.dist]:
                        bin2tips[child.dist].remove(child)
                        if not bin2tips[child.dist]:
                            del bin2tips[child.dist]
                        bin2tips[parent.dist + child.dist].append(child)
                    grandparent.add_child(child, dist=parent.dist + child.dist)
                    grandparent.remove_child(parent)
        # if while pruning we pruned all the root's children by one, we need to update the top part of the tree
        while len(tree.children) == 1:
            tree = tree.children[0]
            tree.dist = 0
            tree.up = None

        # save the pruned tree
        tree.write(outfile=output_tree)
        # save the sequence ids of the pruned tree
        with open(output_ids, 'w+') as f:
            f.write('\n'.join(_.name for _ in tree))


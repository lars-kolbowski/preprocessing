import pandas as pd
import numpy as np
import networkx as nx
import time
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from joblib import Parallel, delayed
import sys

#package modules
import AveragineModel as AM
import ProteoFileReader as PFR

#constants
H_MASS = 1.007276466583
CA_DISTANCE = 1.0033548378


def mz_to_df(mz, intensity, min_idx=0, max_idx=-1):
    """
    Given MZ and Intensity returns a dataframe.
    """
    df = pd.DataFrame()
    df["mz"] = mz
    df["intensity"] = intensity
    if max_idx == -1:
        max_idx = df.shape[0]
    df_filt = df.iloc[min_idx:max_idx]
    df_filt.to_clipboard(sep="\t", index=False, header=False)
    return(df_filt)

def plot_G(G):
    """

    :param G: graph
    :return:
    """
    pos = hierarchy_pos(G, sorted(list(G.nodes()))[0])
    nx.draw(G, pos=pos, with_labels=True)
    edge_labels = dict([((u, v, ), d['label']) for u, v, d in G.edges(data=True)])
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.show()


def hierarchy_pos(G, root):
    '''If there is a cycle that is reachable from root, then result will not be a hierarchy.

       G: the graph
       root: the root node of current branch
       width: horizontal space allocated for this branch - avoids overlap with other branches
       vert_gap: gap between levels of hierarchy
       vert_loc: vertical location of root
       xcenter: horizontal location of root

    Copy & Paste from stackoverflow. Thanks folks!
    '''

    def h_recur(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5,
                pos=None, parent=None, parsed=[]):
        if(root not in parsed):
            parsed.append(root)
            if pos == None:
                pos = {root:(xcenter, vert_loc)}
            else:
                pos[root] = (xcenter, vert_loc)
            neighbors = G.neighbors(root)

            if len(neighbors) != 0:
                dx = width/len(neighbors)
                nextx = xcenter - width/2 - dx/2
                for neighbor in neighbors:
                    nextx += dx
                    pos = h_recur(G, neighbor, width=dx, vert_gap=vert_gap,
                                  vert_loc=vert_loc-vert_gap,
                                  xcenter=nextx, pos=pos, parent=root,
                                  parsed=parsed)
        return pos
    return h_recur(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5)

def plot_graph(G):
    """
    Plot the graph with a random graph layout
    """
    pos = nx.random_layout(G)
    nx.draw(G, pos)
    edge_labels = dict([((u, v, ), d['label']) for u, v, d in
                         G.edges(data=True)])
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    nx.draw_networkx_labels(G, pos, labels={i:i for i in G.nodes()},
                            font_size=16)
    plt.show()


class Deisotoper():
    """
    Description.
    """
    def __init__(self, ppm_tolerance=20, min_score=0.6, min_charge=1, max_charge=7,
                 min_abundance=0.25, min_improve=0.3, verbose=False):
        """
        Initializes the Deisotoper class.

        :param ppm: Maximum deviation in ppm error from an expected m/z to the measured m/z.
        :param allow_ambiguous: Allows the addition of ambiguous charge states.
        :param min_score:  Minimum score for keeping isotope clusters.
        :param min_charge: Int, minimum charge
        :param max_charge: Int, maximum charge stte

        #TODO:
        #think about a useful class & parameter handling.
        
        
        Usage:
        -------------------------------
        Usage is very simple but limited to MGFs at the moment.
        
            >>> import Deisotoper as DE
            >>> infile = "data/test.mgf"
            >>> dt = DE.Deisotoper()
            >>> deisotoped_spectra = dt.deisotope_spectra(infile, 
                                                          show_progress=True,
                                                          n_jobs=-1)
        """
        self.ppm_tolerance = ppm_tolerance
        self.min_score = min_score
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.min_abundance = min_abundance
        self.min_improve = min_improve
        self.verbose = verbose

    
    def parallel_helper(self, spectrum, return_type, show_progress=False,
                        ndone=-1):
            """
            Work horse of the class.
            
            Parameters:
            -------------------------
            return_type: string,
                         either "df" or "spectrum". Df will return a pandas
                         dataframe and spectrum a full featured mgf
                         with charge annotation.
            """
            if ndone % 5000 == 0:
                print("{} spectra done".format(ndone))
                
            mz = spectrum.getMZ()
            intensity = spectrum.getIntensities()
            
            #create graph
            try:
                G = self.spec2graph(mz, intensity, self.ppm_tolerance,
                                    self.min_charge, self.max_charge)
                
                #no clusters in the spectrum continue
                if len(G) == 0:
                    return(spectrum)
                    
                #extract all path with possible isotope clusters
                cluster_ar = self.extract_isotope_cluster(G, mz, verbose=self.verbose)
                
                #resolve ambiguous cluster with scoring
                cluster_res = self.resolve_ambiguous(cluster_ar, mz, intensity,
                                                     min_score=self.min_score,
                                                     min_abundance=self.min_abundance,
                                                     min_improve=self.min_improve,
                                                     verbose=self.verbose)
                
                #write results to dataframe
                cluster_df = self.assemble_spectrum(cluster_res, mz, intensity, 
                                                    spectrum.getTitle()) 
            except:
                print (spectrum.getTitle())
                idwrite = np.random.randint(0, 1000)
                print("Writing erroneous file: {}".format("Error_MGF_{}".format(idwrite)))
                print(spectrum)
                with open("Error_MGF_{}".format(idwrite), 'w') as fobj:
                    fobj.write(spectrum.to_mgf())
                sys.exit()
                
                
            if return_type == "df":
                return(cluster_df)
            else:
                spectrum.peaks = np.array([cluster_df["mz"].values, 
                                           cluster_df["intensity"].values]).transpose()
                spectrum.peakcharge = cluster_df["charge"]
                return(spectrum)
                
                
    def deisotope_spectra(self, infile, in_type="MGF", n_jobs=-1, 
                          return_type="spectrum", show_progress=False):
        """
        Function to deisotope spectra
        """
                
        #process a MGF
        if in_type=="MGF":
            MGF_file = PFR.MGF_Reader()
            MGF_file.load(infile)
            
            results_store = Parallel(n_jobs=n_jobs)\
                (delayed(self.parallel_helper)(spectrum, return_type, show_progress, ii) for ii, spectrum in enumerate(MGF_file))
        
        else:
            sys.exit("Other Types not yet supported. Please provide a MGF")
            
        if return_type == "df":
            results_store_df = pd.concat(results_store)
            return(results_store_df)
            
        else:
            return(results_store)
        
       
    @staticmethod
    def spec2graph(mz, intensity, ppm_error=20, min_charge=1, max_charge=6):
        """
        Creates a graph representation of the peakdata. Connect all peaks
        that are only within a 'da_error' distance measurement.

        :param max_charge:
        :param min_charge:
        :param mz:
        :param intensity:
        :param ppm_error:
        :return:
        """
        #init graph
        n = len(mz)
        G = nx.DiGraph()
       
        zrange = np.arange(float(min_charge), float(max_charge), 1)
#        distance_charge_map = {round((CA_DISTANCE/i), 4): i for i in zrange}
#        charge_distances_map = {j: i for i, j in distance_charge_map.items()}
#        zdist_values = np.array(list(charge_distances_map.values()))

        #build a graph from the spectrum
        #add edges if the distances between peaks matches any isotope distance for charges
        #from min_charge to max_charge
        ambiguity_dic = {}
        for i in range(0, n):
            for j in range(i+1, n):
                # mz difference between next and current peak
                expected_mz = mz[i] + CA_DISTANCE / zrange
                mz_diff = mz[j] - expected_mz
                mz_error_ppm = np.abs((mz_diff / expected_mz)*10**6)

                #only one charge state can match the next peak
                min_idx = mz_error_ppm.argmin()
                if mz_error_ppm[min_idx] <= ppm_error:
                    
                    edgei = (i, j,
                                  {"intensity": intensity[j],
                                   "label": "z:{} \n dist:{} \n ppm:{} \n mz:{}".format(
                                               min_idx+1, np.round(mz_diff[min_idx],4),
                                               np.round(mz_error_ppm[min_idx],2),
                                               np.round(mz[i],2)), "charge": zrange[min_idx]})                    
                    #helper, store only the edge with the lowest mass error
                    #for the same charge
                    if i in ambiguity_dic:
                        #check which entry is better
                        #now check if the current edge is better than the last
                        if ambiguity_dic[i]["error"][min_idx] <= np.round(mz_error_ppm[min_idx],2):
                            pass
                        else:
                            ambiguity_dic[i]["edges"][min_idx] = edgei
                    else:
                        #init new entry
                        ambiguity_dic[i] = {"error":np.zeros_like(zrange),
                                             "edges": [None]*len(zrange)}
                        
                    ambiguity_dic[i]["error"][min_idx] = np.round(mz_error_ppm[min_idx],2)
                    ambiguity_dic[i]["edges"][min_idx]  = edgei
                    
        #only write the best edges to the graph
        edges = []
        for entry in ambiguity_dic.keys():
            for edge in ambiguity_dic[entry]["edges"]:
                if edge is not None:
                    edges.append(edge)
        G.add_edges_from(edges)
        return (G)

    def extract_isotope_cluster(self, G, mz, verbose=False):
        """

        :param G: graph
        :param mz: mz array
        :return:
            
            
        mz_to_df(mz, intensity, min_idx=0, max_idx=-1)
        
        """
        #%%

        #return obj
        isotope_clusters = []

        # start by extracting the ambiguous clusters (connected components)
        # each component can have multiple charge states
        # ideally, we can resolve them by averagine.
        cluster_count = 0
        for subgraph in nx.weakly_connected_component_subgraphs(G):
            has_path = True
            #store the paths here
            path_dic = {}
            while has_path:
                subgraph_copy = subgraph.copy()
                
                start_node = sorted(list(subgraph_copy.nodes()))[0]
                if verbose:
                    print ("Start Graph: ", subgraph.nodes())
                    print ("Checking start node: ", start_node)

                #init node is the first one in the graph
                for start_edge in subgraph_copy[start_node]:

                    #current path defined by charge
                    current_z = subgraph_copy[start_node][start_edge]["charge"]
                    has_charge_path = True
                    current_node = start_edge
                    path = [start_node, start_edge]
                    #while you find an edge with charge = current_z -> follow the path

                    while has_charge_path:
                        if verbose:
                            plot_graph(subgraph)
                            print("Node:", current_node)
                            print("    Graph: ", sorted(subgraph.nodes()))
                            print("    MZ: ", mz[sorted(subgraph.nodes())])
                            print("    Edges: ", subgraph.edges())
                            print("    Path: ", path)
                            print("    CurrentNode: ", current_node)
                            print("    New Charge:", current_z)

                        #plot_graph(subgraph)
                        #print ("Nodes:", subgraph.nodes())
                        next_node = [node for node in subgraph_copy[current_node] if subgraph_copy[current_node][node]["charge"] == current_z]
                        #branch = [node for node in subgraph_copy[current_node] if subgraph_copy[current_node][node]["charge"] != current_z]
                        branch = [edgei for edgei in subgraph.out_edges(current_node, data=True)  if edgei[2]["charge"] != current_z] + \
                                [edgei for edgei in subgraph.in_edges(current_node, data=True)  if edgei[2]["charge"] != current_z]
                        if len(next_node) == 0:
                            has_charge_path = False


                        else:
                            path.append(next_node[0])
                            if len(branch) == 0:
                                subgraph.remove_node(current_node)
                            else:
                                subgraph.remove_edge(current_node, next_node[0])
                                
                            if verbose:
                                print("    REM Node: {}".format(next_node[0]))
                                print("    REM Edge: {}-{}".format(current_node, next_node[0]))

                            current_node = next_node[0]
                    path_dic[(current_z, cluster_count)] = path
                    cluster_count += 1

                subgraph.remove_node(start_node)
                if verbose:
                    print("    Removing Start Node: {} ".format(start_node))
                    print("Final Path: ", path_dic)
                    print("#######################################")

                if len(subgraph) == 0:
                    has_path = False
            isotope_clusters.append(path_dic)
        #%%
        return(isotope_clusters)

    @staticmethod
    def resolve_ambiguous(cluster_ar, mz, intensity, min_score=0.6,
                          min_abundance=0.25, min_improve=0.3,
                          verbose=False, GT=None):
        """

        Function that resolves ambiguous peak clusters. The function uses
        sevaral heuristics that are not optimal. Description of the decision
        process:

        case 1: There is a single connected component (isotope cluster (IC))
                -> do the averagine calculation and score the measured
                intensity vs. the expected intensity. if the score is higher
                than min_score keep the cluster.

        case 2: There are multiple IC where only 1 clusters contributes
                more than 10p of the intensity (Lasso fitting)
                --> Ignore all clusters with <10p abundance and keep the
                left over cluster if case 1 is fullfilled.

        case 3: There are multiple IC where  2 or more clusters contribute
                more than 10p of the intensity (Lasso fitting).
                --> Score the isotopes before intensity correction. Then
                use the fitted Lasso coefficients to substract the overlapping
                peaks and extract "clean" isotope clusters.

                case 3.1: The fitting improved all the scores.
                --> Keep IC if case 1 fullfilled

                case 3.2: The fitting

        Parameters:
        --------------
        cluster_ar:   ar-like,
                    array of isotope clusters that were created by
                    extract_isotope_clusters. <(z,count)>:<mz_peakids>.
        mz: ar-like,
            numpy array of mzs
        intensity: ar-like,
                numpy array of intensities
        min_score: float,
                    minimum pearson score for a averagine fit to be included.
        min_abundance: foat,
                     constraint for resolving ambiguous clusters. The
                     abundance is estimated via Lasso regression.
        min_improve: float
                     Summed delta of all isotope scores needs to be better
                     than this score. E.g. If after the intensity correction
                     of overlapping isotope clusters the score is not better
                     by min_improve reject the correction and score the initial
                     clusters.
        verbose: bool,
                 verbose printing enabled?
        GT: dictionary,
            a Groundtruth dictionary used for debugging and annotating plots.

                    .
        """
        #%%
        result_dic = {}
        nc = 0
        single_cluster = 0
        for cluster_dic in cluster_ar:
            nc += 1
            score_dic = {}
            score_opt_dic = {}
            if verbose:
                #print ("{} / {} Done.".format(nc, n))
                print ("_________________________________________________________________________________________________________________")
                print ("Resolving cluster:", cluster_dic)
                print ("MZ Range: ", np.round(mz[list(cluster_dic.values())[0]][0],2))
            cluster_ids = list(cluster_dic.keys())
            #single cluster does not need any correction, just score it and leave it be
            if len(cluster_dic) == 1:
                single_cluster += 1
                for idx, (z,count) in enumerate(cluster_ids):

                    score = AM.score_cluster(mz[cluster_dic[(z,count)]], intensity[cluster_dic[(z,count)]], z)[0]
                    score_dic[(z,count)] = score
                    result_dic[(z, count)] = (True if score >= min_score else \
                                             False, cluster_dic[(z,count)])


                    if verbose:
                        print ("\t Single Cluster {}. Score: {}".format("Rejected" if score < min_score else "Accepted",
                               np.round(score, 2)))
                continue

            #more than one clusters
            else:

                #here the interesting part starts, resolving overlapping clusters
                # number of paths defines the dimensions
                dims = len(cluster_ids)
                #mapping from nodes to indices
                all_idx = np.sort(np.array(list({item for vals in cluster_dic.values() for item in vals})))
                idx_mapping = {i:j for i,j in zip(all_idx, np.arange(0, len(all_idx), 1))}
                #initialize X matrix for Lasso, columns will correspond to the response factor
                #based on the averagine model
                res = np.zeros(dims*len(all_idx)).reshape(dims, len(all_idx))


                # =============================================================
                #  Step 1: prescoring without optimization
                # =============================================================
                for idx, z in enumerate(cluster_ids):
                    score, peaks = AM.score_cluster(mz[cluster_dic[z]], intensity[cluster_dic[z]], z[0])
                    score_dic[z] = score
                    res[idx][[idx_mapping[i] for i in cluster_dic[z]]] = peaks[:, 1]

                # =============================================================
                # Step 2: Perform non-negative Lasso
                # =============================================================
                #perform, non-negative lasso
                X = np.transpose(res)
                y = intensity[list(all_idx)] #0.0001
                lin = Lasso(alpha=1,precompute=True,max_iter=1000, positive=True, random_state=9999,
                            selection='random', fit_intercept=False,
                            tol=0.001).fit(X,y)

                coefs = (lin.coef_ + 0.000001 )
                abundance_estimate = coefs / coefs.sum()

                # =============================================================
                #  Step 3: Correct Intensities
                # =============================================================
                #coefs will give an estimate of the abundance of the individual
                #species
                abundant_clusters = np.where(abundance_estimate >= min_abundance)[0]
                if len(abundant_clusters) == 1:

                    #do not do any correction, probably just an artefact
                    # e.g. overlapping clusters from z=2 and z=1
                    # 500, 500.5, 501, 501.5 and 500, 501
                    id_tmp = cluster_dic[cluster_ids[abundant_clusters[0]]]
                    score, _ = AM.score_cluster(mz[id_tmp], intensity[id_tmp],
                                                cluster_ids[abundant_clusters[0]][0])
                    score_dic[cluster_ids[abundant_clusters[0]]] = score
                    result_dic[cluster_ids[abundant_clusters[0]]] = (True \
                    if score >= min_score else False, id_tmp)

                    if verbose:
                        print ("\t Multi-Cluster (LA) {}. Score: {}. Cluster: {} Abundance: {}".format(
                                "Rejected" if score < min_score else "Accepted",
                                np.round(score, 2),
                                cluster_ids[abundant_clusters[0]],
                                abundance_estimate))
                    continue

                #okay, here we have two abundant fragments
                #do the correction
                #rescore with optimized intensities
                if verbose:
                    exp_ar = []
                    mz_ar = []
                    z_ar = []


                for idx, z in enumerate(cluster_ids):
                    idx_c = [idx_mapping[idi] for idi in cluster_dic[z]]
                    X_tmp = np.copy(X)
                    X_tmp[:, [i for i in np.arange(0, len(cluster_dic)) if i == idx]] = 0.0

                    #the true intensity of this cluster is the
                    #observed intensity - predicted intensity of the other peps
                    exp_intensity = np.abs(intensity[all_idx] - lin.predict(X_tmp))[idx_c]
                    score, _ = AM.score_cluster(mz[cluster_dic[z]], exp_intensity, z[0])
                    score_opt_dic[z] = score
                    if verbose:
                        print ("\t Multi-Cluster (HA):")
                        print ("\t\t",";".join([str(i) for i in
                                         np.round(abundance_estimate, 2)]))
                        print ("\t\t ({}, {}) {} Score: {}".format(idx, z,
                               "Rejected" if score < min_score else "Accepted",
                               score))

                    if verbose:
                        exp_ar.append(exp_intensity)
                        mz_ar.append(mz[cluster_dic[z]])
                        z_ar.append(z)

            # =================================================================
            # Step 4: Final selection of clusters
            # =================================================================
            #%%
            old_scores = np.array([score_dic[i] for i in score_dic.keys()])
            new_scores = np.array([score_opt_dic[i] for i in score_dic.keys()])
            diff = (new_scores - old_scores).sum()

            if diff < min_improve:
                if verbose:
                    print("Did't improve...")
                #lasso and modelling didn't improve the clusters
                # check
                for i in score_dic.keys():
                    result_dic[i] = (True if score_dic[i] >= min_score \
                                    else False, cluster_dic[i])

            else:
                if verbose:
                    print("Improved by lasso...")
                    print(old_scores)
                    print(new_scores)
                #jay, we improved the scores 'significantly' (#not)
                #lets keep all the clusters we found
                for i in score_opt_dic.keys():
                    result_dic[i] = (True if score_opt_dic[i] >= min_score \
                                    else False, cluster_dic[i])

            if GT:
                print()
                print("Optimized: ", score_dic)
                print("Unoptimized:", score_opt_dic)
                #avoid the ratio to contain a zero
                ratio = [str(i) for i in np.round(coefs / coefs.min(),2)]
                offset = 0.15
                plt.bar(mz[all_idx], intensity[all_idx], label="measured",
                        width=.15, hatch="//")
                for mzi, inti, zi in zip(mz_ar, exp_ar, z_ar):
                    plt.bar(mzi+offset, inti, width=0.15, label="Z: {} Score: {} Score opt: {}".format(zi[0],
                            np.round(score_dic[zi],2),
                            np.round(score_opt_dic[zi], 2)))
                    offset += 0.15
                plt.legend()


                plt.title("""Expected ratio: {} \n
                              Computed ratio: {} \n
                              Case: {}""".format(GT["ratio"], ":".join(ratio),
                                                 GT["Case"]))

                plt.savefig("testcase_{}_resolved.png".format(str.zfill(GT["TestID"], 2)))
            #%%
        return(result_dic)


    @staticmethod
    def assemble_spectrum(cluster, mz, intensity, spectrum_id):
        """
        Creates a dataframe with charge and cluster annotation.
        
        Parameter:
        ---------------
        cluster: array,
                    of clusters
        mz: array
            mz array of a single spectrum
        intensity: array,
                intensity array of a single spectrum
        spectrum: <any>,
                 can be a string or number referencing the spectrum id
        """
        #%%
        ambiguous = []
        charge_states = np.zeros_like(mz)
        cluster_ids = np.zeros_like(mz)
        cluster_count = 0
        for key, (passed_cluster, mz_idx) in cluster.items():
            #if the cluster passed the scoring
            if passed_cluster:
                cluster_count += 1
                
                # check if the cluster is ambiguous
                for peak_id in mz_idx:
                    # if zero, then it is stil free
                    if cluster_ids[peak_id] == 0:
                        cluster_ids[peak_id] = cluster_count
                        charge_states[peak_id] = key[0]
                    else:
                        #not free any longer, need another field
                        ambiguous.append([peak_id, cluster_count, key[0], 
                                          mz[peak_id], intensity[peak_id]])
    
        if len(ambiguous) == 0:
            spectrum_df = pd.DataFrame()
            spectrum_df["mz"] = mz 
            spectrum_df["intensity"] = intensity
            spectrum_df["charge"] = charge_states
            spectrum_df["cluster_id"] = cluster_ids
            spectrum_df["SpectrumID"] = spectrum_id
        else:
            peaks_new, cluster_id_new, charge_new, mz_new, int_new = zip(*ambiguous)
            spectrum_df = pd.DataFrame()
            spectrum_df["mz"] = np.append(mz, np.array(mz_new))
            spectrum_df["intensity"] = np.append(intensity, np.array(int_new))
            spectrum_df["charge"] = np.append(charge_states, np.array(charge_new))
            spectrum_df["cluster_id"] = np.append(cluster_ids, np.array(cluster_id_new))
            spectrum_df["SpectrumID"] = spectrum_id
        #%%
        return(spectrum_df)



"""
Max p regionalization

Heuristically form the maximum number (p) of regions given a set of n
areas and a floor constraint.
"""

__author__ = "Serge Rey <srey@asu.edu>, David Folch <dfolch@fsu.edu>"


#import pysal
from components import check_contiguity
import copy
import numpy as np
from pysal.region import randomregion as RR
import time
from multiprocessing import Pool, Manager
import multiprocessing
from functools import partial

__all__ = ["Maxp", "Maxp_LISA"]

LARGE = 10 ** 6
MAX_ATTEMPTS = 100
lck = multiprocessing.Lock()


class Maxp:
 
    def __init__(self, w, z, floor, floor_variable, verbose=False, initial=100, seeds=[]):
        self.w = w
        self.z = z
        self.floor = floor
        self.floor_variable = floor_variable
        self.verbose = verbose
        self.seeds = seeds
        #file = open("/home/vsindhu1/runtimep.txt","w")
        start_time_in = time.time()
        self.regions, self.area2region, self.p = self.initial_solution()
        #file.write("--- %s seconds for initial solution---" % (time.time() - start_time_in)) 
        print("--- %s seconds for initial solution---" % (time.time() - start_time_in))
        if not self.p:
            self.feasible = False
        else:
            self.feasible = True
            best_val = self.objective_function()
            self.current_regions = copy.copy(self.regions)
            self.current_area2region = copy.copy(self.area2region)
            #self.current_borders = copy.copy(self.borders)
            #self.current_neighbors_b = copy.copy(self.neighbors_b)
            self.initial_wss = []
            self.attempts = 0
            #print("Initial Best Value: ",best_val)
            start_time_loop = time.time()
            #print("Viney")
            p = Pool(64)
            results = p.map(self.th_func,range(initial))
            p.close()
            p.join()
            #print(results)
            #print("B Best Value: ",type(best_val))
            #print("B Regions: ",self.current_regions)
            #print("B A2R: ",self.current_area2region)
#            alist = []
#            blist = []
#            for i in range(initial):
#                alist.append(results[i][3])
#                blist.append(len(results[i][0]))
#            aindex = alist.index(sorted(alist)[0])
#            bestsol = results[aindex]
#            print('Best solution: ',alist)
#            print('Number of reg: ',blist)
#            self.current_regions = bestsol[0]
#            self.current_area2region = bestsol[1]
#            best_val = bestsol[3]
#            for i in range(initial):
#                self.regions = results[i][0]
#                self.area2region = results[i][1]
#                self.p = len(self.regions)
#                best_val = results[i][3]
#                print('Best Value: ',best_val)
#                print('Number of regions: ',self.p)
#                self.swap()
#                print('New Best Value: ',self.objective_function())
#            print('Best Values: ',alist)
#            print("Reg: ",c_regions)
#            print("A2R1: ",c_a2r)
            for i in range(initial):
                if results[i][2] >= self.p:
                    if results[i][3] < best_val:
                        self.current_regions = results[i][0]
                        self.current_area2region = results[i][1]
                        self.p = results[i][2]
                        best_val = results[i][3]
            print("Best Value: ",best_val)
#            print("Regions: ",self.current_regions)
#            print("A2R: ",self.current_area2region)
            #results = self.th_func()
            #print(type(results[0][5]))
            #print(results)
            #print("Viney")
            #file.write("--- %s seconds for loop---" % (time.time() - start_time_loop))
            print("--- %s seconds for loop---" % (time.time() - start_time_loop))
            self.regions = copy.copy(self.current_regions)
            self.p = len(self.regions)
            self.area2region = self.current_area2region
            #self.borders = copy.copy(self.current_borders)
            #self.neighbors_b = copy.copy(self.current_neighbors_b)
            print("Number of regions", self.p)
#            if verbose:
#                print("smallest region ifs: ", min([len(region) for region in self.regions]))
#                print("---Viney Mulle---")
#
#                raw_input='wait'
#				
            start_time_swap = time.time()

            self.swap()
            #file.write("--- %s seconds for swap---" % (time.time() - start_time_swap))
            print("--- %s seconds for swap---" % (time.time() - start_time_swap),"New Value: ", self.objective_function())
            #file.close()

    def th_func(self,i):
        #lck.acquire()
        regions, a2r, p = self.initial_solution()
        if self.initial_solution():
            val = self.objective_function(regions)
        return regions, a2r, p, val
        #lck.release()

    def initial_solution(self):
        self.p = 0
        solving = True
        attempts = 0
        while solving and attempts <= MAX_ATTEMPTS:
            regions = []
            borders = []
            neighbors_b = []
            enclaves = []
            if not self.seeds:
                process = multiprocessing.current_process()
                candidates = copy.copy(self.w.id_order)
                candidates = np.random.RandomState(seed=int(process.pid+time.time())).permutation(candidates)
                candidates = candidates.tolist()
#                print('Candidates: ',process.pid,candidates)
            else:
                seeds = copy.copy(self.seeds)
                nonseeds = [i for i in self.w.id_order if i not in seeds]
                candidates = seeds
                candidates.extend(nonseeds)
            while candidates:
                st= time.time()
                seed = candidates.pop(0)
                # try to grow it till threshold constraint is satisfied
                region = [seed]
                building_region = True
                while building_region:
                    # check if floor is satisfied
                    if self.check_floor(region):
                        regions.append(region)
                        building_region = False
                    else:
                        potential = []
                        for area in region:
                            neighbors = self.w.neighbors[area]
                            neighbors = [neigh for neigh in neighbors if neigh in candidates]
                            neighbors = [neigh for neigh in neighbors if neigh not in region]
                            neighbors = [neigh for neigh in neighbors if neigh not in potential]
                            potential.extend(neighbors)
                        if potential:
                            # add a random neighbor
                            neigID = np.random.randint(0, len(potential))
                            neigAdd = potential.pop(neigID)
                            region.append(neigAdd)
                            # remove it from candidates
                            candidates.remove(neigAdd)
                        else:
                            #print 'enclave'
                            #print region
                            enclaves.extend(region)
                            building_region = False
                print('Each process time: ', process.pid, time.time() - st)
            # check to see if any regions were made before going to enclave stage
            if regions:
                feasible = True
            else:
                attempts += 1
                break
            self.enclaves = enclaves[:]
            a2r = {}
            for r, region in enumerate(regions):
                for area in region:
                    a2r[area] = r
            encCount = len(enclaves)
            encAttempts = 0
            while enclaves and encAttempts != encCount:
                enclave = enclaves.pop(0)
                neighbors = self.w.neighbors[enclave]
                neighbors = [neighbor for neighbor in neighbors if neighbor not in enclaves]
                candidates = []
                for neighbor in neighbors:
                    region = a2r[neighbor]
                    if region not in candidates:
                        candidates.append(region)
                if candidates:
                    # add enclave to random region
                    regID = np.random.randint(0, len(candidates))
                    rid = candidates[regID]
                    regions[rid].append(enclave)
                    a2r[enclave] = rid
                    # structure to loop over enclaves until no more joining is possible
                    encCount = len(enclaves)
                    encAttempts = 0
                    feasible = True
                else:
                    # put back on que, no contiguous regions yet
                    enclaves.append(enclave)
                    encAttempts += 1
                    feasible = False
            #for reg in regions:
                #neig = []
                #for area in reg:
                  #  neigh = self.w.neighbors[area]
                 #   neig.append(neigh)
                    #print("Area",area)
                    #print("Neigh",neigh)
                #cnt = 0
                #bordr = []
                #print("Neig: ",neig)
               # neig11 = []
                
                #for l in neig:
                  #  ll = [x for x in l if x not in reg]
                 #   neig11.append(ll)
                #sys.stdout.write('Neighbors: %s \n' % (neig11))
                #for area in reg:
                    #if neig11[cnt]:
                    #    #print("Neig Count",neig[cnt])
                   #     #print("Count",cnt)
                  #      bordr.append(area)
                 #   cnt+=1
                #neig1 = sum(neig11,[])
                #neig2 = list(set(neig1))
                #print("Neig2",neig2)
                #print("Bordr",bordr)
                #neighbors_b.append(neig2)
                #borders.append(bordr)
            #print("Border Areas: ",borders)
            #print("Neighbor Areas: ",neighbors_b)
            if feasible:
                solving = False
                #self.regions = regions
                #self.area2region = a2r
                #self.p = len(regions)
                p = len(regions)
                #self.borders = borders
                #self.neighbors_b = neighbors_b
                return regions, a2r, p
            else:
                if attempts == MAX_ATTEMPTS:
                    print('No initial solution found')
                    self.p = 0
                attempts += 1
                return False

    def swap(self):
        swapping = True
        swap_iteration = 0
        if self.verbose:
            print('Initial solution, objective function: ', self.objective_function())
            sys.stdout.write('Regions: %s \n' % (self.regions))
            sys.stdout.write('Best Val: %f' % (self.best_val))
        total_moves = 0
        self.k = len(self.regions)
        #print("Total regions:", self.k)
        changed_regions = [1] * self.k
        nr = list(range(self.k))
        #print('Regions: %s \n' % (self.regions))
        manager = Manager()
        shared_list = manager.list()
        shared_lock = manager.list()
        smoves = manager.Value('i',0)
#        print('Initial moves: ',smoves.value)
        for i in range(self.k):
            shared_lock.append(manager.Lock())
        for region in self.regions:
            shared_list.append(region)
        shared_map = manager.list()
        for reg in self.area2region:
            shared_map.append(reg)
        #print("Shared List: ", shared_list)
        #shared_list[0].append(13)
        #print("Shared1 List: ", shared_list)
#        file = open("/home/vsindhu1/runtimeswap.txt","w")
        
        while swapping:
            self.moves_made = 0
            regionIds = [r for r in nr if changed_regions[r]]
            np.random.permutation(regionIds)
            #print("Region Ids: ",regionIds)
            #changed_regions = [0] * self.k
            swap_iteration += 1
            pst = time.time()
            p = Pool(self.k)
#            print('Process Start time: ',time.time() - pst)
            mpsttm = time.time()
            results_swap = p.map(partial(self.swapfunc, shared_list=shared_list,shared_lock=shared_lock,smoves=smoves), range(self.k))
#            print("Total iteration time: ", time.time() - mpsttm)
            p.close()
            p.join()
            #self.regions = []
            #print("Result Swap: ", results_swap)
            for i in range(self.k):
                self.moves_made+=results_swap[i]
            #self.regions = results_swap[0][1]
            #print("Updated Shared List: ", shared_list)
            self.regions = shared_list
            for r, region in enumerate(self.regions):
                for area in region:
                    self.area2region[area] = r
            #print('objective function: ', self.objective_function())
            #print('moves_made: ', self.moves_made)
            total_moves += self.moves_made
            #print('Total moves_made: ', total_moves)
            #print('After Swap Regions: %s \n' % (self.regions))
#            print('objective function: ', self.objective_function())
            if self.moves_made == 0:
                swapping = False
                self.swap_iterations = swap_iteration
                self.total_moves = total_moves
            #if self.verbose:
                print('Number of iterations: ', self.swap_iterations)
        print('Total Moves: ',smoves.value)
                #print('objective function: ', self.objective_function())
        
        

    def swapfunc(self,seed,shared_list,shared_lock,smoves):
        #shared_list[seed].append(self.regions[seed])
#        file1 = open("/home/vsindhu1/runtimeswap.txt","a")
        psttm = time.time()
        members = shared_list[seed]
        m_made = 0
        neighbors = []
        for member in members:
            candidates = self.w.neighbors[member]
            candidates = [candidate for candidate in candidates if candidate not in members]
            candidates = [candidate for candidate in candidates if candidate not in neighbors]
            neighbors.extend(candidates)
        candidates = []
        try:
            for neighbor in neighbors:
                block = copy.copy(shared_list[self.area2region[
                    neighbor]])
                if check_contiguity(self.w, block, neighbor):
                    block.remove(neighbor)
                    fv = self.check_floor(block)
                    if fv:
                        candidates.append(neighbor)
        except:
            pass
        # find the best local move
        if candidates:
            nc = len(candidates)
            moves = np.zeros([nc, 1], float)
            best = None
            cv = 0.0
            try:
                for area in candidates:
                    current_internal = shared_list[seed]
                    current_outter = shared_list[self.area2region[
                        area]]
                    current = self.objective_function([current_internal, current_outter])
                    new_internal = copy.copy(current_internal)
                    new_outter = copy.copy(current_outter)
                    new_internal.append(area)
                    new_outter.remove(area)
                    new = self.objective_function([new_internal,
                                                   new_outter])
                    change = new - current
                    if change < cv:
                        best = area
                        cv = change
            except:
                pass
            if best:
                # make the move
                lcksttm = time.time()
                old_region = self.area2region[best]
                lock1 = shared_lock[old_region]
                lock2 = shared_lock[seed]
                if old_region < seed:
                    lock1.acquire()
                    lock2.acquire()
                else:
                    lock2.acquire()
                    lock1.acquire()
#                lck.acquire()
                #print("Started")
                #print("Best: ", best)
                area = best
#                old_region = self.area2region[area]
#                shared_lock[old_region].acquire()
                #print("Area: ", area)
                #print("Old Region: ", old_region)
                #print("A Regions before swap: ", shared_list)
                try:
                    ls = copy.copy(shared_list[self.area2region[
                    neighbor]])
                    fv1 = self.check_floor(ls)
                    if fv1:
                        so_list = shared_list[old_region]
                        so_list.remove(area)
                        shared_list[old_region] = so_list
                        self.area2region[area] = seed
                        ss_list = shared_list[seed]
                        ss_list.append(area)
                        shared_list[seed] = ss_list
                        m_made += 1
                        smoves.value += 1
                except:
                    pass
                #print("A Regions after swap: ", shared_list)
                #changed_regions[seed] = 1
                #changed_regions[old_region] = 1
                #print("Ended")
#                lck.release()
#                shared_lock[old_region].release()
#                shared_lock[seed].release()
#                shared_lock[old_region].acquire()
                if old_region < seed:
                    lock2.release()
                    lock1.release()
                else:
                    lock1.release()
                    lock2.release()
#                file1.write("Total Lock Time: %s\n" % str(seed))
#                file1.write("%s \n" % (time.time() - lcksttm))
#                print("Total Lock Time: ",seed,": ", time.time() - lcksttm)
#        print("Total Process Time: ",seed,": ", time.time() - psttm)
#        print("Number of areas: ",seed,": ",len(shared_list[seed]))
#        file1.write("Total Process Time: %s\n" % str(seed))
#        file1.write("%s \n" % (time.time() - psttm))
#        file1.close()
        return m_made


    def check_floor(self, region):
        selectionIDs = [self.w.id_order.index(i) for i in region]
        cv = sum(self.floor_variable[selectionIDs])
        if cv >= self.floor:
            #print len(selectionIDs)
            return True
        else:
            return False

    def objective_function(self, solution=None):
        # solution is a list of lists of region ids [[1,7,2],[0,4,3],...] such
        # that the first region has areas 1,7,2 the second region 0,4,3 and so
        # on. solution does not have to be exhaustive
        if not solution:
            solution = self.regions
        wss = 0
        for region in solution:
            selectionIDs = [self.w.id_order.index(i) for i in region]
            m = self.z[selectionIDs, :]
            var = m.var(axis=0)
            wss += sum(np.transpose(var)) * len(region)
        return wss

    def inference(self, nperm=99):
        """Compare the within sum of squares for the solution against
        simulated solutions where areas are randomly assigned to regions that
        maintain the cardinality of the original solution.

        Parameters
        ----------

        nperm       : int
                      number of random permutations for calculation of
                      pseudo-p_values

        Attributes
        ----------

        pvalue      : float
                      pseudo p_value

        Examples
        --------

        Setup is the same as shown above except using a 5x5 community.

        >>> import numpy as np
        >>> import pysal
        >>> np.random.seed(100)
        >>> w=pysal.weights.lat2W(5,5)
        >>> z=np.random.random_sample((w.n,2))
        >>> p=np.ones((w.n,1),float)
        >>> floor=3
        >>> solution=pysal.region.Maxp(w,z,floor,floor_variable=p,initial=100)

        Set nperm to 9 meaning that 9 random regions are computed and used for
        the computation of a pseudo-p-value for the actual Max-p solution. In
        empirical work this would typically be set much higher, e.g. 999 or
        9999.

        >>> solution.inference(nperm=9)
        >>> solution.pvalue
        0.1

        """
        ids = self.w.id_order
        num_regions = len(self.regions)
        wsss = np.zeros(nperm + 1)
        self.wss = self.objective_function()
        cards = [len(i) for i in self.regions]
        sim_solutions = RR.Random_Regions(ids, num_regions,
                                          cardinality=cards, permutations=nperm)
        cv = 1
        c = 1
        for solution in sim_solutions.solutions_feas:
            wss = self.objective_function(solution.regions)
            wsss[c] = wss
            if wss <= self.wss:
                cv += 1
            c += 1
        self.pvalue = cv / (1. + len(sim_solutions.solutions_feas))
        self.wss_perm = wsss
        self.wss_perm[0] = self.wss

    def cinference(self, nperm=99, maxiter=1000):
        """Compare the within sum of squares for the solution against
        conditional simulated solutions where areas are randomly assigned to
        regions that maintain the cardinality of the original solution and
        respect contiguity relationships.

        Parameters
        ----------

        nperm       : int
                      number of random permutations for calculation of
                      pseudo-p_values

        maxiter     : int
                      maximum number of attempts to find each permutation

        Attributes
        ----------

        pvalue      : float
                      pseudo p_value

        feas_sols   : int
                      number of feasible solutions found

        Notes
        -----

        it is possible for the number of feasible solutions (feas_sols) to be
        less than the number of permutations requested (nperm); an exception
        is raised if this occurs.

        Examples
        --------

        Setup is the same as shown above except using a 5x5 community.

        >>> import numpy as np
        >>> import pysal
        >>> np.random.seed(100)
        >>> w=pysal.weights.lat2W(5,5)
        >>> z=np.random.random_sample((w.n,2))
        >>> p=np.ones((w.n,1),float)
        >>> floor=3
        >>> solution=pysal.region.Maxp(w,z,floor,floor_variable=p,initial=100)

        Set nperm to 9 meaning that 9 random regions are computed and used for
        the computation of a pseudo-p-value for the actual Max-p solution. In
        empirical work this would typically be set much higher, e.g. 999 or
        9999.

        >>> solution.cinference(nperm=9, maxiter=100)
        >>> solution.cpvalue
        0.1

        """
        ids = self.w.id_order
        num_regions = len(self.regions)
        wsss = np.zeros(nperm + 1)
        self.cwss = self.objective_function()
        cards = [len(i) for i in self.regions]
        sim_solutions = RR.Random_Regions(ids, num_regions,
                                          cardinality=cards, contiguity=self.w,
                                          maxiter=maxiter, permutations=nperm)
        self.cfeas_sols = len(sim_solutions.solutions_feas)
        if self.cfeas_sols < nperm:
            raise Exception('not enough feasible solutions found')
        cv = 1
        c = 1
        for solution in sim_solutions.solutions_feas:
            wss = self.objective_function(solution.regions)
            wsss[c] = wss
            if wss <= self.cwss:
                cv += 1
            c += 1
        self.cpvalue = cv / (1. + self.cfeas_sols)
        self.cwss_perm = wsss
        self.cwss_perm[0] = self.cwss


class Maxp_LISA(Maxp):
    """Max-p regionalization using LISA seeds

    Parameters
    ----------

    w              : W
                     spatial weights object
    z              : array
                     nxk array of n observations on k variables used to
                     measure similarity between areas within the regions.
    y              : array
                     nx1 array used to calculate the LISA statistics and
                     to set the intial seed order
    floor          : float
                     value that each region must obtain on floor_variable
    floor_variable : array
                     nx1 array of values for regional floor threshold
    initial        : int
                     number of initial feasible solutions to generate
                     prior to swapping

    Attributes
    ----------

    area2region     : dict
                      mapping of areas to region. key is area id, value is
                      region id
    regions         : list
                      list of lists of regions (each list has the ids of areas
                      in that region)
    swap_iterations : int
                      number of swap iterations
    total_moves     : int
                      number of moves into internal regions


    Notes
    -----

    We sort the observations based on the value of the LISAs. This
    ordering then gives the priority for seeds forming the p regions. The
    initial priority seeds are not guaranteed to be separated in the final
    solution.

    Examples
    --------

    Setup imports and set seeds for random number generators to insure the
    results are identical for each run.

    >>> import numpy as np
    >>> import pysal
    >>> np.random.seed(100)

    Setup a spatial weights matrix describing the connectivity of a square
    community with 100 areas.  Generate two random data attributes for each area
    in the community (a 100x2 array) called z. p is the data vector used to
    compute the floor for a region, and floor is the floor value; in this case
    p is simply a vector of ones and the floor is set to three. This means
    that each region will contain at least three areas.  In other cases the
    floor may be computed based on a minimum population count for example.

    >>> w=pysal.lat2W(10,10)
    >>> z=np.random.random_sample((w.n,2))
    >>> p=np.ones(w.n)
    >>> mpl=pysal.region.Maxp_LISA(w,z,p,floor=3,floor_variable=p)
    >>> mpl.p
    30
    >>> mpl.regions[0]
    [99, 98, 89]

    """
    def __init__(self, w, z, y, floor, floor_variable, initial=100):

        lis = pysal.Moran_Local(y, w)
        ids = np.argsort(lis.Is)
        ids = ids[list(range(w.n - 1, -1, -1))]
        ids = ids.tolist()
        mp = Maxp.__init__(
            self, w, z, floor=floor, floor_variable=floor_variable,
            initial=initial, seeds=ids)


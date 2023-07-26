#include<bits/stdc++.h>
#include <cstdlib>
#include <ctime>

using namespace std;

#define Max_Capacity 200

int Num_Customers;
int Num_Chromosomes;
int alpha;
int beta;

class Customer {
public:
    int cust_id;
    float x_coord, y_coord, demand, ready_time, due_time, service_time;


    Customer(int cust_id, float x_coord, float y_coord, float demand, float ready_time, float due_time, float service_time) {
        this->cust_id = cust_id;
        this->x_coord = x_coord;
        this->y_coord = y_coord;
        this->demand = demand;
        this->ready_time = ready_time;
        this->due_time = due_time;
        this->service_time = service_time;

    }
    void print_data() {
        cout << this->cust_id << " " << this->x_coord << " " << this->y_coord << " " << this->demand << " " << this->ready_time << " " << this->due_time << " " << this->service_time << endl ;
    }
};

vector<Customer> customers;

class Cost {

public:
    int num_vehicles;
    float dist;
    float fitness;
    vector<int>routes;


    Cost(int num_vehicles, float dist, float fitness, vector<int>routes) {
        this->num_vehicles = num_vehicles;
        this->dist = dist;
        this->fitness = fitness;
        this->routes = routes;
    }

};

void take_input() {
    freopen("inputs/input_c103.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    for (int i = 0; i < Num_Customers; i++) {
        int cust_id;
        float x_coord, y_coord, demand, ready_time, due_time, service_time;
        cin >> cust_id >> x_coord >> y_coord >> demand >> ready_time >> due_time >> service_time;
        cust_id -= 2;
        Customer cust(cust_id, x_coord, y_coord, demand, ready_time, due_time, service_time);
        customers.push_back(cust);
    }
}

// returns a vector of random numbers of length len in range of 'range'
vector<int> random_vec(int range, int len) {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    int r;
    unordered_set<int> st;
    while (st.size() < len) {
        r = std::rand();
        r = r % range;
        st.insert(r);
    }
    vector<int> randoms;
    for (auto x : st) {
        randoms.push_back(x);
        //cout<<x<<endl;
    }
    return randoms;
}


vector<vector<int>> Generate_initial_population() {
    vector<int> init(Num_Customers);
    for (int i = 0; i < Num_Customers; i++) init[i] = i;
    set<vector<int>> st;
    st.insert(init);
    while (st.size() < Num_Chromosomes) {
        random_shuffle(init.begin(), init.end());
        st.insert(init);
    }
    vector<vector<int>> init_pop;
    for (auto v : st)
    {
        init_pop.push_back(v);
    }
    random_shuffle(init_pop.begin(), init_pop.end());
    return init_pop;
}


vector<int> pareto_ranking(vector<pair<int, float>> costs) {
    vector<int> rank(costs.size(), INT_MAX);
    int r = 1;
    while (1) {
        bool all_done = true;
        for (int i = 0; i < costs.size(); i++) {
            if (rank[i] == INT_MAX) {
                all_done = false;
                bool flag = true;
                // checking if i is non-dominated
                for (int j = 0; j < costs.size(); j++) {
                    if (rank[j] == INT_MAX) {
                        if (costs[j].first <= costs[i].first and costs[j].second <= costs[i].second) {
                            if (costs[j].first < costs[i].first or costs[j].second < costs[i].second) {
                                flag = false;
                                break;
                            }
                        }
                    }
                }
                if (flag) rank[i] = r;
            }

        }
        if (all_done) break;
        r++;
    }

    return rank;

}

Cost cost_function(vector<int> chromosome) {
    int X = 40, Y = 50;
    float prev_x = X, prev_y = Y;
    float time = 0, dist = 0, demand = 0;
    vector<int>routes = {0};
    // starting a route from pos 0
    //routes.push_back(0);
    int num_vehicles = 1;
    for (int i = 0; i < chromosome.size(); i++) {

        Customer gene = customers[chromosome[i]];
        float a = gene.x_coord;
        float b = gene.y_coord;
        float d = sqrt((prev_x - a) * (prev_x - a) + (prev_y - b) * (prev_y - b));
        if (demand + gene.demand <= Max_Capacity and (time + d) <= gene.due_time) { // removed + service time from (time)
            time += d;
            if (time < gene.ready_time) {
                time = gene.ready_time;
            }
            float x = gene.x_coord;
            float y = gene.y_coord;

            float dist_bw = sqrt((prev_x - x) * (prev_x - x) + (prev_y - y) * (prev_y - y));
            dist += dist_bw;

            demand += gene.demand;

            time += gene.service_time;

            prev_x = x;
            prev_y = y;
        }
        else {

            // going back to warehouse
            // starting a new route from pos i
            routes.push_back(i);
            num_vehicles++;
            time = 0;
            demand = 0;
            prev_x = X;
            prev_y = Y;

            // adding dist of last customer in previous route to origin
            i--;

            gene = customers[chromosome[i]];
            float x = gene.x_coord;
            float y = gene.y_coord;

            float dist_bw = sqrt((prev_x - x) * (prev_x - x) + (prev_y - y) * (prev_y - y));
            dist += dist_bw;
        }
    }

    // adding dist of last customer to origin
    prev_x = X;
    prev_y = Y;
    Customer gene = customers[chromosome[chromosome.size() - 1]];
    float x = gene.x_coord;
    float y = gene.y_coord;

    float dist_bw = sqrt((prev_x - x) * (prev_x - x) + (prev_y - y) * (prev_y - y));
    dist += dist_bw;

    float fitness = alpha * num_vehicles + beta * (dist);
    Cost cost(num_vehicles, dist, fitness, routes);

    return cost;
}
// <------- MUTATIONS ------->

vector<int> Inversion_mutation(vector<int> chromosome, int rand_index) {

    auto child = chromosome;

    //swap(child[rand_index], child[rand_index + 1]);
    reverse(child.begin() + rand_index, child.begin() + rand_index + 3);
    return child;
}

vector<int> Swap_mutation(vector<int> chromosome, int rand_index1, int rand_index2) {
    auto child = chromosome;
    swap(child[rand_index1], child[rand_index2]);
    return child;
}
vector<int> Cycle_mutation(vector<int> chromosome, vector<int>random_k) {
    //replace gene with cycle rotation
    int temp = chromosome[random_k[0]];
    for (int i = 1; i < random_k.size(); i++) {
        chromosome[random_k[i - 1]] = chromosome[random_k[i]];
    }
    chromosome[random_k[random_k.size() - 1]] = temp;
    return chromosome;
}
void Mutation(vector<vector<int>> &population) {



    int pop_10 = population.size() / 10;
    vector<int> rand_pop = random_vec(population.size() - 1, pop_10);
    vector<int> rand_indexes = random_vec(population[0].size() - 2, pop_10);
    vector<int> rand_indexes2 = random_vec(population[0].size(), pop_10); //for swap mutation
    vector<int> random_k = random_vec(population[0].size() - 1, 4);


    for (int i = 0; i < pop_10; i++) {
        //swap
        //auto child = Swap_mutation(population[rand_pop[i]], rand_indexes[i], rand_indexes2[i]);

        // cyclic mutation
        //auto child = Cycle_mutation(population[rand_pop[i]], random_k);

        // inversion mutation
        auto child = Inversion_mutation(population[rand_pop[i]], rand_indexes[i]);


        population.push_back(child);
    }


}

// <----- CROSSOVERS ------>

// BCRC CROSSOVER
vector<int> bcrc_util(vector<int> &del_par, unordered_set<int>& stp, int nv, float fit, float dist) {
    //(stp.begin(), stp.end());

    vector<int> new_par;
    for (auto gene : stp) {
        vector<int> temp;
        float temp_fit = fit;
        int temp_nv = nv;
        float temp_dist = dist;
        for (int i = 0; i <= del_par.size(); i++) {
            // for (int j = 0; j < i; j++) v.push_back(del_par[j]);
            // v.push_back(gene);
            // for (int j = i; j < del_par.size(); j++) v.push_back(del_par[j]);
            vector<int> v = del_par;
            v.insert(v.begin() + i, gene);

            Cost c_v = cost_function(v);
            if (c_v.dist <= temp_dist and c_v.fitness <= temp_fit) {
                temp = v;
                temp_fit = c_v.fitness;
                temp_dist = c_v.dist;
                temp_nv = c_v.num_vehicles;
            }
        }
        if (temp.size() == 0) {
            del_par.push_back(gene);
        }
        else del_par = temp;
    }
    //cout << del_par.size() << endl;
    return del_par;
}
vector<vector<int>> bcrc(vector<vector<int>>& population, int fath, int moth) {

    Cost cf = cost_function(population[fath]);
    Cost cm = cost_function(population[moth]);

// random indexes in routes vector
    int rf = random_vec(cf.routes.size() - 1, 1)[0];
    int rm = random_vec(cm.routes.size() - 1, 1)[0];

    int rfs = cf.routes[rf];
    int rms = cm.routes[rm];

    int rfe, rme;
    if (rf == cf.routes.size() - 1) {
        rfe = Num_Customers;
    }
    else {
        rfe = cf.routes[rf + 1];
    }

    if (rm == cm.routes.size() - 1) {
        rme = Num_Customers;
    }
    else {
        rme = cm.routes[rm + 1];
    }

    unordered_set<int> stf, stm;
    for (int i = rfs; i < rfe; i++) {
        stf.insert(population[fath][i]);
    }
    for (int i = rms; i < rme; i++) {
        stm.insert(population[moth][i]);
    }

    vector<int> del_fath, del_moth;

    for (int i = 0; i < Num_Customers; i++) {
        if (stm.find(population[fath][i]) == stm.end()) {
            del_fath.push_back(population[fath][i]);
        }
    }

    for (int i = 0; i < Num_Customers; i++) {
        if (stf.find(population[moth][i]) == stf.end()) {
            del_moth.push_back(population[moth][i]);
        }
    }

    auto ch1 = bcrc_util(del_fath, stm, cf.num_vehicles, cf.fitness, cf.dist);
    auto ch2 = bcrc_util(del_moth, stf, cm.num_vehicles, cm.fitness, cm.dist);


    return {ch1, ch2};
}


// PMX CROSSOVER

vector<int> pmx_util(vector<int> p1, vector<int> p2, int x, int y) {
    int range = p1.size();
    vector<int> ch(range, -1);
    unordered_map<int, int> mp1, mp2;
    unordered_set<int> st;
    for (int i = 0; i < range; i++) {
        mp2[p2[i]] = i;
        mp1[p2[i]] = p1[i];
    }

    for (int i = x; i <= y; i++) {
        ch[i] = p1[i];
        st.insert(p1[i]);
    }


    for (int i = x; i <= y; i++) {
        if (st.count(p2[i])) continue;
        int a = p2[i], b = mp1[p2[i]];

        while (1) {
            int pos = mp2[b];
            //cout<<"a-> "<<a<<"  b-> "<<b<<endl;
            if (ch[pos] != -1) {
                b = ch[pos];
                //cout<<pos<<" "<<ch[pos]<<endl;
            }
            else {
                ch[pos] = a;
                break;
            }
        }
    }
    for (int i = 0; i < range; i++) {
        if (ch[i] == -1) ch[i] = p2[i];
    }

    return ch;
}
vector<vector<int>> pmx(vector<int> p1, vector<int> p2) {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    int range = p1.size();
    int x, y;
    x = (std:: rand()) % range;
    y = (std:: rand()) % range;
    while (x == y) {
        y = (std:: rand()) % range;
    }
    if (x > y) swap(x, y);
    vector<int> ch1 = pmx_util(p1, p2, x, y);
    vector<int> ch2 = pmx_util(p2, p1, x, y);
    return {ch1, ch2};
}

// Cyclic Crossover
void cyclic_util(int idx, vector<int>& par1, vector<int>& par2, vector<int> &ch, unordered_map<int, int> &mp) {
    if (ch[idx] != -1) return ;
    ch[idx] = par1[idx];
    cyclic_util(mp[par2[idx]], par1, par2, ch, mp);
}
vector<vector<int>> cyclic(vector<int> &par1, vector<int> &par2) {
    int len = par1.size();
    vector<int> ch1(len, -1), ch2(len, -1);
    unordered_map<int, int> mp1, mp2;
    for (int i = 0; i < len; i++) {
        mp1[par1[i]] = i;
        mp2[par2[i]] = i;
    }
    bool flag = true;
    for (int i = 0; i < len; i++) {
        if (ch1[i] != -1) continue;
        if (flag) {
            cyclic_util(i, par1, par2, ch1, mp1);
        }
        else {
            cyclic_util(i, par2, par1, ch1, mp2);
        }
        flag = !flag;
    }
    flag = false;
    for (int i = 0; i < len; i++) {
        if (ch2[i] != -1) continue;
        if (flag) {
            cyclic_util(i, par1, par2, ch2, mp1);
        }
        else {
            cyclic_util(i, par2, par1, ch2, mp2);
        }
        flag = !flag;
    }
    return {ch1, ch2};
}
void Crossover(vector<vector<int>>& population) {
    int x = (population.size() * 8) / 10;
    //int x = population.size();
    if (x % 2 == 1) x--;
    vector<int> random_parents = random_vec(population.size() - 1, x);
    //cout << random_parents.size() << " random parents" << endl;

    for (int i = 0; i < random_parents.size(); i = i + 2)
    {
        //cout << i << endl;
        //auto children = BCRC_Crossover(population, random_parents[i], random_parents[i + 1]);
        //auto children = pmx(population[random_parents[i]], population[random_parents[i]]);
        //auto children = cyclic(population[random_parents[i]], population[random_parents[i]]);
        //cout << children.size() << endl;
        auto children = bcrc(population, random_parents[i], random_parents[i + 1]);
        population.push_back(children[0]);
        population.push_back(children[1]);
    }
}




int main() {

    Num_Customers = 100;
    Num_Chromosomes = 100;

    alpha = 500;
    beta = 1;

    take_input();

    vector<vector<int>> init_pop = Generate_initial_population();

    int gen = 0;
    while (gen < 300) {
        auto temp = init_pop;
        Crossover(temp);
        Mutation(temp);

        vector<pair<int, float>> costs(temp.size());
        for (int i = 0; i < temp.size(); i++) {
            auto  cost = cost_function(temp[i]);
            costs[i] = {cost.num_vehicles, cost.dist};
            // cout << cost.num_vehicles << " " << cost.dist << " " << cost.fitness << endl;
        }

        auto rank = pareto_ranking(costs);

        vector<pair<int, vector<int>>> vpp(temp.size());
        for (int i = 0; i < temp.size(); i++) {
            vpp[i] = {rank[i], temp[i]};
        }
        sort(vpp.begin(), vpp.end());

        int percent_worst = 0;
        for (int i = 0; i < Num_Customers - percent_worst; i++) {
            init_pop[i] = vpp[i].second;
        }

        for (int i = Num_Customers - percent_worst; i < Num_Customers; i++) {
            init_pop[i] = vpp[vpp.size() - 1 - percent_worst - i + Num_Customers].second;
        }

        gen++;
    }


    vector<pair<int, float>> costs(init_pop.size());
    for (int i = 0; i < init_pop.size(); i++) {
        auto  cost = cost_function(init_pop[i]);
        costs[i] = {cost.num_vehicles, cost.dist};
        // cout << cost.num_vehicles << " " << cost.dist << " " << cost.fitness << endl;
    }

    auto rank = pareto_ranking(costs);
    int best_num_veh = 100, avg_num_veh = 0;
    float best_dist = 10000, avg_dist = 0;
    vector<int> best_dist_chrom;
    vector<int> best_vehicle_chrom;
    int cnt = 0;
    for (int i = 0; i < init_pop.size(); i++) {
        if (rank[i] == 1) {
            //cout << costs[i].first << " " << costs[i].second << endl;
            cnt++;
            best_num_veh = min(best_num_veh, costs[i].first);
            avg_num_veh += costs[i].first;
            //best_dist = min(best_dist, costs[i].second);
            if (best_dist > costs[i].second) {
                best_dist = costs[i].second;
                best_dist_chrom = init_pop[i];
            }
            avg_dist += costs[i].second;

        }
    }
    avg_num_veh /= cnt;
    avg_dist /= cnt;
    cout << "Best Num-Vehicles -> " << best_num_veh << "\n";
    cout << "Avg Num-Vehicles -> " << avg_num_veh << "\n";
    cout << "Best Distance -> " << best_dist << "\n";
    cout << "Avg Distance -> " << avg_dist << "\n";
    cout << "chromosome with best dist -> " << "\n";
    // for (int i = 0; i < best_dist_chrom.size(); i++) {
    //     cout << best_dist_chrom[i] << " ";
    // }
    // cout << "\n";
    auto cv = cost_function(best_dist_chrom);
    cout << "Num of vehicle -> " << cv.num_vehicles << " Distance -> " << cv.dist << endl;
    int i = 0;
    int j = 1;
    cv.routes.push_back(100);
    while (i < 100) {
        cout << "Route " << j << " ->  ";
        while (i < 100 and i < cv.routes[j]) {
            cout << best_dist_chrom[i++] << " ";
        }
        j++;
        cout << endl;
    }




    //vector<int> v = {2, 21, 73, 41, 56, 4, 5, 83, 61, 85, 37, 93, 14, 44, 38, 43, 13, 27, 69, 76, 79, 3, 54, 24, 80, 28, 12, 40, 53, 26, 30, 51, 9, 66, 1, 31, 88, 7, 10, 33, 29, 78, 34, 35, 77, 36, 47, 19, 8, 46, 17, 39, 23, 67, 55, 25, 45, 82, 18, 84, 60, 89, 52, 6, 59, 99, 94, 96, 62, 11, 90, 20, 32, 70, 63, 64, 49, 48, 65, 71, 81, 50, 68, 72, 75, 22, 74, 58, 92, 42, 15, 87, 57, 97, 95, 98, 16, 86, 91, 100};
    //vector<int> v = {83, 45, 61, 84, 5, 60, 89, 87, 57, 2, 58, 28, 29, 78, 34, 35, 3, 77, 73, 22, 72, 54, 24, 80, 12, 42, 15 , 41, 75, 56, 74, 21, 94, 96, 99, 6, 65, 71, 81, 20, 32, 70, 50, 33, 30, 51, 9, 66, 1, 40, 53, 92, 37, 98 , 85, 16, 86, 91, 97 , 13, 62, 88, 8, 46, 17, 93, 59, 63, 64, 90, 10, 31 , 27, 69, 76, 79, 68, 52, 7, 11, 19, 49, 48, 82, 14, 44, 38, 43, 100, 95, 36, 47, 18, 39, 23, 67, 55, 4, 25, 26};
    //vector<int> v = { 92 , 98 , 14 , 44 , 38 , 86 , 16 , 61 , 85 , 91 , 100 , 37, 42 , 43 , 15 , 57 , 41 , 74 , 72 , 73 , 21 , 58, 7 , 19 , 11 , 8 , 46 , 47 , 48 , 82 , 18 , 89, 71 , 65 , 78 , 34 , 35 , 81 , 77 , 28, 26 , 39 , 23 , 67 , 55 , 24 , 29 , 3, 27 , 69 , 30 , 9 , 66 , 20 , 51 , 1, 2 , 22 , 75 , 56 , 4 , 25 , 54, 36 , 64 , 49 , 63 , 90 , 32 , 70, 52 , 62 , 88 , 84 , 17 , 93 , 59 , 50 , 33 , 76 , 79 , 10 , 31, 60 , 45 , 83 , 5 , 99 , 6 , 94 , 96 , 95 , 97 , 87 , 13 , 40 , 53 , 12 , 68 , 80};

    // for (int i = 0; i < 100; i++) v[i]--;
    // auto cv = cost_function(v);
    // cout << cv.num_vehicles << " " << cv.dist << endl;

    // int i = 0;
    // int j = 1;
    // cv.routes.push_back(100);
    // while (i < 100) {
    //     while (i < 100 and i < cv.routes[j]) {
    //         cout << v[i++] << " ";
    //     }
    //     j++;
    //     cout << endl;
    // }


    // for (auto v : init_pop) {
    //     for (auto x : v) {
    //         cout << x << " ";
    //     }
    //     cout << endl;
    // }


}
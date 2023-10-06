#ifndef MATERSDK_LEVEL_H
#define MATERSDK_LEVEL_H


#include <vector>
#include <utility>
#include <algorithm>


namespace matersdk {
namespace mtp {


class Combinations {
public:
    Combinations(std::vector<std::vector<std::pair<int, int>>> mjus_njus_lst, bool sort_unique_mark=false);

    void show() const;

    const std::vector<std::vector<std::pair<int, int>>> get_combinations() const;

    const int get_num_combinations() const;

    static int get_level(const int mju, const int nju);

private:
    std::vector<std::vector<std::pair<int, int>>> mjus_njus_lst;
};  // class: Combinations




class CombinationsSortBasis {
public:
    CombinationsSortBasis(const Combinations rhs) : combinations(rhs)
    {}


    bool operator()(int index_i, int index_j) const {
        // Step 1. 计算 index_i 的 mtp level
        int level_index_i = 0;
        for (int ii=0; ii<combinations.get_combinations()[index_i].size(); ii++) {
            level_index_i += Combinations::get_level(
                                    combinations.get_combinations()[index_i][ii].first, 
                                    combinations.get_combinations()[index_i][ii].second);
        }

        // Step 2. 计算 index_j 的 mtp level
        int level_index_j = 0;
        for (int ii=0; ii<combinations.get_combinations()[index_j].size(); ii++) {
            level_index_j += Combinations::get_level(
                                    combinations.get_combinations()[index_j][ii].first,
                                    combinations.get_combinations()[index_j][ii].second);
        }

        return (level_index_i < level_index_j);
    }

private:
    Combinations combinations;
};  // class : CombinationsSortBasis



class CombinationsArrangement {
public:
    CombinationsArrangement(const Combinations combinations, int* new_indices) : combinations(combinations), new_indices(new_indices)
    {}

    Combinations arrange() const {
        std::vector<std::vector<std::pair<int, int>>> sorted_mjus_njus_lst;
        sorted_mjus_njus_lst.clear();
        for (int ii=0; ii<this->combinations.get_combinations().size(); ii++) {
            sorted_mjus_njus_lst.push_back(combinations.get_combinations()[this->new_indices[ii]]);
        }

        Combinations sorted_combinations(sorted_mjus_njus_lst);
        return sorted_combinations;
    };

private:
    Combinations combinations;
    int* new_indices;
};



class MTPLevel {
public:
    MTPLevel(int mju, int nju);

    // Combinations: mju, nju
    void calc_combinations(int num_M, int max_level);

private:
    int mju;
    int nju;
};  // class: MTPLevel







Combinations::Combinations(std::vector<std::vector<std::pair<int, int>>> mjus_njus_lst, bool sort_unique_mark) {
    // Step 1. Init the Combinations with `std::vector<std::vector<std::pair<int, int>>>`
    this->mjus_njus_lst.resize(mjus_njus_lst.size());
    for (int ii=0; ii<mjus_njus_lst.size(); ii++) 
        this->mjus_njus_lst[ii].resize(mjus_njus_lst[ii].size());
    
    for (int ii=0; ii<mjus_njus_lst.size(); ii++) {
        for (int jj=0; jj<mjus_njus_lst[ii].size(); jj++) {
            this->mjus_njus_lst[ii][jj] = mjus_njus_lst[ii][jj];
        }
    }

    // Step 2. 
    if (sort_unique_mark == true) {

    }

    // Step 3. 
}


void Combinations::show() const {
    for (int ii=0; ii<this->mjus_njus_lst.size(); ii++) {
        printf("[mju, nju] :\t");
        int level = 0;
        for (int jj=0; jj<this->mjus_njus_lst[ii].size(); jj++) {
            printf(
                "[%4d, %4d],  ", 
                this->mjus_njus_lst[ii][jj].first, 
                this->mjus_njus_lst[ii][jj].second
            );
            level += Combinations::get_level(this->mjus_njus_lst[ii][jj].first, this->mjus_njus_lst[ii][jj].second);
        }
        printf("level = %4d,\tno.%3d\n", level, ii);
    }
}


const std::vector<std::vector<std::pair<int, int>>> Combinations::get_combinations() const {
    return (const std::vector<std::vector<std::pair<int, int>>>)this->mjus_njus_lst;
}


const int Combinations::get_num_combinations() const {
    return this->mjus_njus_lst.size();
}


int Combinations::get_level(const int mju, const int nju) {
    return (2 + 4*mju + nju);
}





MTPLevel::MTPLevel(int mju, int nju) {
    this->mju = mju;
    this->nju = nju;
}


void MTPLevel::calc_combinations(int num_M, int max_level) {

}


}; // namespace : mtp
}; // namespace : matersdk
#endif 
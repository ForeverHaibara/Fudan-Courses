#include<iostream>
#include<algorithm>
#include<vector>
#include<random>
#include<time.h>
using namespace std;

//将 vec[left: right) 中小于 x 的数放在数组左侧，大于等于 x 的数放在数组右侧
//返回小于 x 的数的个数
//注：C++ STL 有现成的 partition函数
template<class T>
inline int partition_(typename vector<T>::iterator left,typename vector<T>::iterator right,T x){
    typename vector<T>::iterator l = left, aim;
    --right;
    while (left < right){
        while (left < right && *left < x){
            ++left;
        }
        while (left < right && *right >= x){
            --right;
        }
        swap(*left,*right);
    }
    return left-l;
}

//返回 vec[left: right) 中从小到大第 k 个数 (1 <= k <= vec.size())
template <class T>
T selectKth(typename vector<T>::iterator left,typename vector<T>::iterator right,int k){
    if (right - left <= 10){
        sort(left,right);
        return *(left+k-1);
    }
    vector<T> medians;
    typename vector<T>::iterator i;
    for (i = left; right-i>=5; i+=5){
        medians.push_back(selectKth<T>(i,i+5,3));
    }
    if (i != right){
        medians.push_back(selectKth<T>(i,right,(1+(int)(right-i))>>1));
    }
    T x = selectKth<T>(medians.begin(),medians.end(),(1+medians.size())>>1);

    int part = partition_(left,right,x);
    if (part >= k) return selectKth<T>(left,left+part,k);
    else if (part < k) return selectKth<T>(left+part,right,k-part);
    else return x; // part == k-1
}

//返回 vec 中从小到大第 k 个数 (1 <= k <= vec.size())
//注：该函数会修改 vec 数组，如果不想修改，请将传引用改为传值
template<class T>
inline T selectKth(vector<T>& vec,int k){
    return selectKth<T>(vec.begin(),vec.end(),k);
}

int main(){
    int n = 10000;
    vector<int> a(n);
    srand(time(NULL));
    for (int i=0;i<n;++i) a[i] = i;
    random_shuffle(a.begin(),a.end());
    cout << "First 10 elements: "; for (int i=0;i<10;++i) cout << a[i] <<' ';cout <<endl;
    cout << selectKth(a,1024) << endl;
    return 0;
}

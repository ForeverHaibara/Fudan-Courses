// convolution with FFT implementation
#include<iostream>
#include<vector>
#include<math.h>
#include<string>
using namespace std;


//********************************************************
//                    Complex Class
//********************************************************
template<class T>
class Complex{
public:
    T re,im;
    Complex(T _r=T(0),T _i=T(0)):re(_r),im(_i){}
    Complex& operator += (const Complex& x){ re += x.re; im += x.im; return *this;}
    Complex& operator -= (const Complex& x){ re -= x.re; im -= x.im; return *this;}
    Complex& operator /= (T        x){ re /= x;    im /= x;    return *this;}
    Complex& operator *= (const Complex& x){
        T _re(re);
        re = re*x.re - im*x.im; 
        im = _re*x.im + im*x.re; 
        return *this;
    }
    operator int(){ return round(re);}
    operator double(){ return (double)re;}
    void Conjugate(){ im = -im; }
    
    template<class T2> friend Complex<T2> operator + (const Complex<T2>& a,const Complex<T2>& b);
    template<class T2> friend Complex<T2> operator - (const Complex<T2>& a,const Complex<T2>& b);
    template<class T2> friend Complex<T2> operator * (const Complex<T2>& a,const Complex<T2>& b);
    template<class T2> friend ostream&    operator <<(ostream& out        ,const Complex<T2>& x);
};

template<class T2>
Complex<T2> operator + (const Complex<T2>& a,const Complex<T2>& b){ 
    return Complex<T2>(a.re+b.re, a.im+b.im); 
}

template<class T2>
Complex<T2> operator - (const Complex<T2>& a,const Complex<T2>& b){ 
    return Complex<T2>(a.re-b.re, a.im-b.im); 
}

template<class T2>
Complex<T2> operator * (const Complex<T2>& a,const Complex<T2>& b){ 
    return Complex<T2>(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re); 
}

template<class T2>
ostream& operator << (ostream& out,const Complex<T2>& x){
    if (x.im >= T2(0)) out << x.re << '+' << x.im << 'i';
    else out << x.re << x.im << 'i';
    return out;
}




//********************************************************
//                Fast Fourier Transform
//********************************************************

// upper2exp(n) return an integer which is the smallest, yet no less than n, exponential of 2
inline int upper2exp(int n){
    int m = 1;
    while (m < n) m <<= 1;
    return m;
}

// fill a vector ``a`` with zeros at the end such that the length is an exponential of 2
template<class T>
inline void fill2exp(vector<T>& a,int n=0){
    if (n < a.size()) n = upper2exp(a.size());
    for (int i=n-a.size()-1;i>=0;--i) a.push_back(T(0));
}

// regard the coefficients of a[start:end:step] as a new polynomial, 
// and calculate the fft result given the unit roots
template<class T2>
vector<Complex<double>> FFTRecursion(const vector<T2>& a,int start,int end,int step,
                                    const vector<Complex<double>>& roots){

    if (start + step >= end) return vector<Complex<double>>({Complex<double>(a[start])});

    vector<Complex<double>> ans((end-start-1)/step + 1);
    vector<Complex<double>> even_coeffs = FFTRecursion(a,start,end-step,step<<1,roots);
    vector<Complex<double>> odd_coeffs  = FFTRecursion(a,start+step,end,step<<1,roots);

    int m = ans.size() >> 1;
    for (int i=0,j=m,k=0; i<m; ++i,++j,k+=step){
        ans[i] = even_coeffs[i] + roots[k]*odd_coeffs[i];
        ans[j] = even_coeffs[i] - roots[k]*odd_coeffs[i];
    }

    return ans;
}

// fast fourier transform of a polynomial
// set the argument inv = true  to perform inverse fast fourier transform
template<class T2>
vector<Complex<double>> FFT(vector<T2> a,bool inv = false){
    // guarantee the length is an exponential of 2
    fill2exp(a);
    int n = a.size();

    vector<Complex<double>> roots(n); // unit roots
    if (inv){
        // inverse fft shall use the conjugated unit roots
        for (int i=n>>1,j=i;i;--i,++j){
            roots[j] = roots[i] = Complex<double>(cos(2*M_PI*i/n) , sin(2*M_PI*i/n));
            roots[i].Conjugate();
        }
    }else{
        for (int i=n>>1,j=i;i;--i,++j){
            roots[j] = roots[i] = Complex<double>(cos(2*M_PI*i/n) , sin(2*M_PI*i/n));
            roots[j].Conjugate();
        }
    }
    roots[0] = Complex<double>(1,0);

    if (!inv) return FFTRecursion(a,0,n,1,roots);

    //else if (inv == true)
    vector<Complex<double>> ans = FFTRecursion(a,0,n,1,roots);
    for (int i=n-1;i>=0;--i){
        ans[i] /= n;
    }
    return ans;
}


// Convolution, or Polynomial Multiplication
template<class T>
vector<T> Convolution(vector<T> a,vector<T> b){
    // n is not less than the degree of polynomial a*b
    int n = upper2exp( a.size() + b.size() );
    fill2exp(a,n);
    fill2exp(b,n);

    vector<Complex<double>> temp1 = FFT(a);
    vector<Complex<double>> temp2 = FFT(b);
    for (int i=n-1;i>=0;--i){
        temp1[i] *= temp2[i];
    }

    vector<T> ans(n);
    temp1 = FFT(temp1,true); // inverse FFT
    for (int i=n-1;i>=0;--i){
        ans[i] = T(temp1[i]);
    }
    return ans;
}


// Long Integers Multiplication
string LongIntegersMultiplication(const string& str_a,const string& str_b){
    vector<int> a(str_a.size());
    vector<int> b(str_b.size());
    for (int i=a.size()-1,j=0;i>=0;--i,++j) a[j] = str_a[i]-'0';
    for (int i=b.size()-1,j=0;i>=0;--i,++j) b[j] = str_b[i]-'0';

    vector<int> ans = Convolution(a,b);

    // prevent Access Overflow Error
    for (int i=0;i<10;++i) ans.push_back(0);

    int n = ans.size();
    for (int i=0;i<n;++i){
        if (ans[i]>9){
            int temp = ans[i]%10;
            ans[i] /= 10;
            for (int k=1;ans[i];++k){
                ans[i+k] += ans[i]%10;
                ans[i] /= 10;
            }
            ans[i] = temp;
        }
    }
    
    string ans2;
    int j = n-1;
    while (ans[j] == 0) --j;  // remove leading zeros
    for (;j>=0;--j) ans2 += (char)(ans[j]+'0');
    return ans2;
}


// Test
int main(){
    /*
    vector<int> a({1,2,3,5,-2,4,3});
    vector<Complex<double>> b = FFT(a);
    for (auto point: b){
        cout << point << ' ';
    }cout << endl;
    
    vector<Complex<double>> c = FFT(b,true);
    for (auto point: c){
        cout << point << ' ';
    }cout << endl;

    vector<int> d({-2,3,1,0,1});
    for (auto v: Convolution(a,d)){
        cout << v << ' ';
    }cout << endl;

    cout << LongIntegersMultiplication("195842401208","5212486512919991754");
    cout << "\nDONE\n";
    */
    string numa,numb;
    cin >> numa >> numb;
    cout << LongIntegersMultiplication(numa,numb);
    return 0;
}

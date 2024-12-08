#include <iostream>
#include <algorithm>
using namespace std;
int k = 4;
float epsilon = pow(10, -6);
double function(double x) {
	return 6 * pow(x, 5);
}
double function2(double x) {
	return pow(x, 1.0 / 16.0) * sqrt(1 + pow(x, 2));
}
string doubleToPow10toN(double value) {
	int exponent = 0;
    char value_str[20];
    if (value != 0) {
        exponent = static_cast<int>(std::floor(std::log10(std::abs(value))));
        value /= std::pow(10, exponent);
    }
    sprintf_s(value_str, "%.2f*10^%+d", value, exponent);
    return string(value_str);
}
double findJ(double h, int n, double a, double b, int functionNumber) {
	double sum1 = 0;
	double sum2 = 0;
	for (int i = 1; i < n; i++) {
		double xi = i * h + a;
		if (functionNumber == 1) { sum1 += function(xi); }
		else { sum1 += function2(xi); }
	}
	for (int i = 0; i < n; i++) {
		double xi = i * h + a;
		if (functionNumber == 1) { sum2 += function(xi + h / 2.0); }
		else { sum2 += function2(xi + h / 2.0); }
	}
	double J = 0;
	if (functionNumber == 1) { J = h / 6.0 * (function(a) + function(b) + 2 * sum1 + 4 * sum2); }
	else { J = h / 6.0 * (function2(a) + function2(b) + 2 * sum1 + 4 * sum2); }
	return J;
}
void print_table(double a, double b, int functionNumber) {
	double prev1_K_delta = 0;
	double prev2_K_delta = 0;
	double prev_delta_runge = 0;
	double IntegralTochnoe = 0;
	double delta_theoretical = 0;
	if (functionNumber == 1) { IntegralTochnoe = 1; }
	else if (functionNumber == 2 && a == 0) { IntegralTochnoe = 1.898840171345417; }
	else { IntegralTochnoe = 0; }
	printf("%-5s %-7s %-15s %-15s %-12s\n", "N", "K_delta", "Delta_tochnoe", "Delta_Runge", "Delta_theoretical");
	printf("%s\n", std::string(65, '-').c_str());
	for (int n = 1; n <= 65536; n *= 2) {
		double h = (b - a) / static_cast<double>(n);
		if (functionNumber == 1) {
			double M4 = (b - a) * max(720 * a, 720 * b);
			double c = M4 / double(2880);
			delta_theoretical = c * pow(h, k);
		}
		else {
			delta_theoretical = 0;
		}
		double Jh = findJ(h, n, a, b, functionNumber);
		double Jh2 = findJ(h / 2.0, n * 2, a, b, functionNumber);
		double Jh4 = findJ(h / 4.0, n * 4, a, b, functionNumber);
		double K_delta = (Jh2 - Jh) / (Jh4 - Jh2);
		double delta_runge = (Jh2 - Jh) / (pow(2, k) - 1);
		if (functionNumber == 1) {
			if (n == 1) {
				prev1_K_delta = K_delta;
				prev_delta_runge = delta_runge;
				printf("%-5d %-7s %-15s %-15s %-12s\n", n, "-", doubleToPow10toN((IntegralTochnoe - Jh)).c_str(), "-", doubleToPow10toN(delta_theoretical).c_str());
			}
			else if (n == 2) {
				prev2_K_delta = K_delta;
				printf("%-5d %-7s %-15s %-15s %-12s\n", n, "-", doubleToPow10toN((IntegralTochnoe - Jh)).c_str(), doubleToPow10toN(prev_delta_runge).c_str(), doubleToPow10toN(delta_theoretical).c_str());
				prev_delta_runge = delta_runge;
			}
			else {
				printf("%-5d %-7.2f %-15s %-15s %-12s\n", n, prev1_K_delta, doubleToPow10toN((IntegralTochnoe - Jh)).c_str(), doubleToPow10toN(prev_delta_runge).c_str(), doubleToPow10toN(delta_theoretical).c_str());
				prev1_K_delta = prev2_K_delta;
				prev2_K_delta = K_delta;
				prev_delta_runge = delta_runge;
			}
		}
		else {
			if (n == 1) {
				prev1_K_delta = K_delta;
				prev_delta_runge = delta_runge;
				printf("%-5d %-7s %-15s %-15s %-12s\n", n, "-", "-", "-", "-");
			}
			else if (n == 2) {
				prev2_K_delta = K_delta;
				printf("%-5d %-7s %-15s %-15s %-12s\n", n, "-", "-", doubleToPow10toN(prev_delta_runge).c_str(), "-");
				prev_delta_runge = delta_runge;
			}
			else {
				printf("%-5d %-7.2f %-15s %-15s %-12s\n", n, prev1_K_delta, "-", doubleToPow10toN(prev_delta_runge).c_str(), "-");
				prev1_K_delta = prev2_K_delta;
				prev2_K_delta = K_delta;
				prev_delta_runge = delta_runge;
			}
		}
	}
	cout << endl;
}
int main() {
	print_table(0, 1, 1);
	print_table(0.0, 1.5, 2);
	print_table(0.001, 1.5, 2);
	return 0;
}
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::string;
using std::vector;

class Student {
private:
	string name;
	int age;
	double score;
public:
	Student(string input_name, int input_age, double input_score) {
		name = input_name;
		age = input_age;
		score = input_score;
	}

	void print_info() {
		std::cout << "Name : " << name <<
		", age : " << age << ", score : " << score
		<< std::endl;
	}
};

void readfile(vector<Student> &students, string &filename) {
	std::ifstream infile(filename);
	if (infile.is_open()) {
		string name;
		int age;
		double score;
		while (infile >> name >> age >> score) {
			Student new_student(name, age, score);
			students.push_back(new_student);
		}
		infile.close();
	} else {
		throw std::invalid_argument("Can't open file to read.");
	}
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cout << "Please input the file name!" << std::endl;
	}
	vector<Student> students;
	string filename = argv[1];
	try {
		readfile(students, filename);
	}
	catch (std::invalid_argument &e) {
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		std::cout << "Something wrong happened!" << std::endl;
		return EXIT_FAILURE;
	}
	for (auto x : students) {
		x.print_info();
	}
	return EXIT_SUCCESS;
}


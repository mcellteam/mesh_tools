#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

int main(int argc,char **argv)
{
  if (argc!=3)
  {
    cout << "\nError. Example usage: sample_data sample_num input_filename\n\n";
    exit(1);
  }
  int sample = atoi(argv[1]);
  if (sample<1)
  {
    cout << "\n\nError. Sample rate ( " << sample << " ) is less than 1.\n\n";
    exit(1);
  }
  char str[2048];
  std::ifstream file(argv[2]);
  if (!file.is_open())
  {
    cout << "\n\nError. Can not open file:" << argv[2] << endl;
    exit(1);
  }
  int i = 1;
  while(!file.eof())
  {
    file.getline(str,2048);
    if (i==sample)
    {
      cout << str << endl;
      i=1;
    }
    else
    {
      i++;
    }
  }
  file.close(); 
}

#ifndef BED_PARSING_H
#define BED_PARSING_H

// This class is useful for reading tab-delimited tables of text.
class Row
{
 public:
  std::string const & operator[](std::size_t index) const
    {
      return m_data[index];
    }
  std::size_t size() const
    {
      return m_data.size();
    }
  void readNextRow(std::istream & str)
  {
    std::string line, cell;

    std::getline(str, line);

    std::stringstream lineStream(line);

    m_data.clear();

    while (std::getline(lineStream, cell, '\t')) {
      // Remove spaces within each tab-delimited cell.
      cell.erase(std::remove(cell.begin(), cell.end(), ' '),
		 cell.end());
      // Delete all occurrences of '^M' (aka '\r').
      cell.erase(std::remove(cell.begin(), cell.end(), '\r'),
		 cell.end());
      m_data.push_back(cell);
    }
  }
 private:
  std::vector<std::string> m_data;
};


#endif

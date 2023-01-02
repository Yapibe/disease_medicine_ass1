# Yair Pickholz Berliner 205435357
import pandas as pd


# define data class
class myData:
    def __init__(self, path) -> None:
        self.data = pd.read_csv(path, low_memory=False)

    def num_players_height(self, x, y):
        """
        return number of players with height between x and y including
        :param x: lower bound of height
        :param y: upper bound of height
        :return: number of players with height between x and y including
        """
        return len(self.data[(self.data['height_cm'] >= x) & (self.data['height_cm'] <= y)])

    def df_birthyear(self, year):
        """
        return a dataframe of players born in year
        :param year: year to filter by
        :return: dataframe of players born in year
        """
        # convert year to string
        year = str(year)
        # return dataframe of players born in year
        df = self.data[self.data['dob'].str.contains(year)]
        # return dataframe with columns short_name, long_name
        df = df[['short_name', 'club_name']]
        return df

    def list_sorted(self, col1, col2, k):
        """
        return a list of tuples (col1, col2) sorted by col2 in descending order
        :param col1: column to choose K highest values from
        :param col2: column to sort by
        :param k: number of players to return
        :return: the list of tuples
        """
        # return a dataframe sorted by col1 in descending order
        df = self.data.sort_values(by=[col1], ascending=False)
        # return top k rows of df
        df = df.head(k)
        # sort by col2 in ascending order
        df = df.sort_values(by=[col2], ascending=True)
        return df['short_name'].tolist()

    def tuples_players_by_year(self, x, y):
        """
        return a list of tuples (year, number of players born in year) sorted by year in ascending order
        :param x: earliest year to include
        :param y: latest year to include
        :return: list of tuples (year, number of players born in year) sorted by year
        """
        tup_list = []
        for i in range(x, y + 1):
            tup_list.append((i, self.df_birthyear(i).shape[0]))
        return tup_list

    def mean_std(self, col, string):
        """
        return the mean and standard deviation of column col
        :param col: column to calculate mean and standard deviation of
        :param string: name of player to filter by
        :return: dataframe with mean and standard deviation of column col for players with name string
        """
        string = string + ' '
        df = self.data[self.data['long_name'].str.startswith(string)]
        return df[col].mean(), df[col].std()

    def max_players(self, col):
        """
        return the maximum value of column col
        :param col: column to calculate maximum value of
        :return: most common value of column col
        """
        # return most common value in column col
        return self.data[col].value_counts().index[0]


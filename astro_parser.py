
""" astro_parse.py

Script that will parse the Astrochymist bibliography of molecules
and dump the year, source, molecule, and observation technique to
a CSV file.

Current issues:
    There are no easy ways to separate the molecular formula
    and name from the way the HTML is written. If David
    ever gets around to fixing the HTML, then we can talk about
    getting this done properly.
"""

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
import requests
import periodictable as pt
import re


def clean_formula_string(formula: str):
    for symbol in ["+", "–", "c-", "trans", "n", "E-", "t-", "Z-", "l-", "-", "≡", "="]:
        formula = formula.replace(symbol, "")
    try:
        form_obj = pt.formula(formula)
    except:
        print(formula)
    return form_obj


def get_atom_dict(formula: str):
    form_obj = clean_formula_string(formula)
    atom_dict = {str(atom): number for atom, number in form_obj.atoms.items()}
    return atom_dict


def formula_multiplicity(formula: str) -> bool:
    """
    Determines if the molecule is open or closed shell from its formula.
    """
    if "+" in formula:
        num_elec = -1
    elif "–" in formula:
        num_elec = 1
    else:
        num_elec = 0
    form_obj = clean_formula_string(formula)
    for atom, number in form_obj.atoms.items():
        num_elec += atom.number * number
    return num_elec % 2 != 0


def clean_html_tags(text: str):
    # this comprises HTML tags that need to be removed. Primarily
    # this is for molecular formula
    for match in ["<sup>", "</sup>", "<sub>", "</sub>", 
                  "<b>", "</b>", "<i>", "</i>", "*", "\n",
                  '''<font face="Helvetica" size="4">'''
                 ]:
        text = text.replace(match, "")
    return text


def final_text_clean(text: str) -> str:
    """
    This is just a function that operates as a pipeline for any
    post cleaning required before adding it to the row data.
    """
    # some molecules are not thoroughly cleaned, and so we remove
    # all the text that comes after <br/>
    if "br" in text:
        text = text.split("<br/>")[0]
    return text


def main():
    url = "http://www.astrochymist.org/astrochymist_ism.html"
    response = requests.get(url)

    soup = BeautifulSoup(response.text, 'lxml')  # Parse the HTML as a string

    table = soup.find_all('table')[2]            # Grab the third table, which contains the molecules
    data = list()

    colormapping = {
        "cyan": "Radio",
        "yellow": "UV/Vis",
        "pink": "IR"
    }
    for index, row in enumerate(table.find_all("tr")[1:]):
        row_data = {"formula": "", "year": 0, 
                    "Radio": False, "IR": False, "UV/Vis": False, "authors": list()}
        # loop over columns
        for col_index, column in enumerate(row.find_all(["th", "td"])):
            # first column is the year; we only grab it if it begins with a digit
            # because there are some weird A+E characters
            if col_index == 0:
                try:
                    year = clean_html_tags(str(column.contents[0].contents[0]))
                    if year.isdigit():
                        row_data["year"] = year
                except ValueError:
                    pass
                # this is the year
            # this grabs the molecular formula; in principle the name is also
            # contained after <br/> but we seldom use it anyway
            elif col_index == 1:
                contents = [str(value) for value in column.contents[0]]
                formula_string = list()
                for substring in contents:
                    if str(substring) == "<br/>":
                        break
                    else:
                        substring = clean_html_tags(substring)
                        formula_string.append(substring)
                formula_string = ''.join(formula_string)
                formula_string = final_text_clean(formula_string)
                if not formula_string[0].isdigit():
                    row_data["formula"] = formula_string
                if column.find("font", {"color": "ORANGE"}):
                    row_data["disputed"] = True
            elif col_index == 2:
                # Loop over the references, and use regex to extract
                # the author names
                authors = list()
                references = ''.join([str(val) for val in column.contents])
                # Regex will pick up capital letter followed by period,
                # space, a capital letter and at least three lower case
                authors.extend(
                    re.findall(r"[A-Z]\.\s[A-Za-z]\w{3,}", references)
                )
                # Find all authors not formatted the same way
                authors.extend(
                    re.findall(r"[A-Z]\.[A-Za-z]\w{3,}", references)
                )
                # Remove whitespace
                authors = [author.replace(" ", "") for author in authors]
                # Make sure all authors are formatted the same way
                authors = [" ".join(author.split(".")) for author in authors]
                row_data["authors"] = authors
            # last column is the detection method; this is inferred from the
            # color used for the HTML
            elif col_index == 3:
                techs = list()
                for color, tech in colormapping.items():
                    if column.find("font", {"color": color}):
                        row_data[tech] = True
        # get the atom numbers too
        atom_dict = get_atom_dict(row_data["formula"])
        row_data["total"] = sum(atom_dict.values())
        row_data.update(atom_dict)
        data.append(row_data)

    # Set up pandas dataframe
    df = pd.DataFrame(data=data)
    # remove blank entries
    df = df.loc[(df["formula"] != "")]
    df["radical?"] = df["formula"].apply(formula_multiplicity)

    # Dump to file
    df.to_csv("astrochymist.csv", index=False)
    return df

def time_analysis(df):
    # Perform analysis of detections of time, as well as
    # the number of discoveries per author over time.
    author_df = pd.DataFrame(df["Authors"].values.tolist(), index=df["Year"])
    unique_authors = list()
    # Create a list of unique authors
    for column in author_df:
        for author in author_df[column].unique():
            if author not in unique_authors:
                unique_authors.append(author)
    # Dictionary for holding stuff
    authorship = list()
    tally = list()
    running_tally = {author: 0 for author in unique_authors}
    for index, row in df.iterrows():
        row_data = {author: False for author in unique_authors}
        for author in row["Authors"]:
            row_data[author] = True
        authorship.append(row_data)
    authorship_df = pd.DataFrame(authorship, index=df["Year"])
    # Collapse the dataframe along years
    authorship_df = authorship_df.groupby(authorship_df.index).sum()
    # Cumulative sum of authors over the years
    tally_df = authorship_df.cumsum(axis=0)
    tally_df.index = authorship_df.index
    # Get the ranking of detections by authors
    rank_df = authorship_df.sum().sort_values(ascending=False) 
    authorship_df.to_csv("author_timeseries.csv")
    tally_df.to_csv("author_cumsum.csv")
    rank_df.to_csv("author_ranking.csv")
    return authorship_df, tally_df, rank_df
    
if __name__ == "__main__":
    df = main()

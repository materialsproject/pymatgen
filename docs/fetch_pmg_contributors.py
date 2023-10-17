"""
Use the GitHub API to fetch the contributors of the pymatgen repository
into a pandas DataFrame, sort it and save it as a CSV file.
"""

# %%
from __future__ import annotations

from datetime import datetime

import pandas as pd
import requests
from tqdm import tqdm

from docs.secrets import GH_TOKEN

__author__ = "Janosh Riebesell"
__date__ = "2023-08-11"


# %%
headers = {"Authorization": f"token {GH_TOKEN}"}  # GH personal access token
base_url = "https://api.github.com/repos/materialsproject/pymatgen"

contributors_url = f"{base_url}/stats/contributors"
pull_requests_url = f"{base_url}/pulls?state=closed&sort=created&direction=asc"

contributors_response = requests.get(contributors_url, headers=headers).json()


# %%
gh_user_data = {}

for contributor in tqdm(contributors_response):
    user_url = f"https://api.github.com/users/{contributor['author']['login']}"
    user_data = requests.get(user_url, headers=headers).json()

    if user_data.get("type") == "Bot":  # Exclude bots
        continue
    gh_user_data[contributor["author"]["login"]] = user_data


# %%
merged_prs_per_contributor = {}
for login in tqdm(gh_user_data):
    merged_pr_response = requests.get(pull_requests_url + f"&creator={login}&state=merged", headers=headers)
    merged_prs_per_contributor[login] = merged_pr_response.json()


# %%
# Process contributor data and PRs
contributor_dict = {}
weeks_with_prs_col = "Weeks with merged PR"
n_contribs_col = "Total contributions"

for contributor in tqdm(contributors_response):
    login = contributor["author"]["login"]
    if login not in merged_prs_per_contributor:
        continue
    pr_list = merged_prs_per_contributor[login]

    weeks_with_merged_pr = sum(1 for week in contributor["weeks"] if week["c"] > 0)

    years_with_merged_pr = {datetime.fromtimestamp(week["w"]).year for week in contributor["weeks"] if week["c"] > 0}

    name = gh_user_data[login]["name"]
    contributor_dict[name] = {
        "GitHub username": login,
        weeks_with_prs_col: weeks_with_merged_pr,
        "Years with merged PR": ", ".join(map(str, years_with_merged_pr)),
        "Number years active": len(years_with_merged_pr),
        "Oldest PR": pr_list[0]["html_url"] if pr_list else None,
        n_contribs_col: contributor["total"],
    }

df_contributors = pd.DataFrame(contributor_dict.values()).sort_values(n_contribs_col, ascending=False)

today = f"{datetime.now():%Y-%m-%d}"
df_contributors.to_csv(f"{today}-top-contributors.csv", index=False)


# %%
# %store df_contributors # cache the data


# restore the data
# %store -r df_contributors

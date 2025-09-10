def get_number_of_pn_tagged_svjs(df, wp):

    j0_pn = df["GoodJetsAK80_pNetJetTaggerScore"]
    j1_pn = df["GoodJetsAK81_pNetJetTaggerScore"]
    j2_pn = df["GoodJetsAK82_pNetJetTaggerScore"]
    j3_pn = df["GoodJetsAK83_pNetJetTaggerScore"]

    svj0 = j0_pn > wp
    svj1 = j1_pn > wp
    svj2 = j2_pn > wp
    svj3 = j3_pn > wp

    svj0 = svj0.astype(int)
    svj1 = svj1.astype(int)
    svj2 = svj2.astype(int)
    svj3 = svj3.astype(int)

    n_svjs = svj0 + svj1 + svj2 + svj3

    return n_svjs


def get_number_of_tagged_svjs(df, wp):

    def __is_tagged(
        wnae_pt_0_200,
        wnae_pt_200_300,
        wnae_pt_300_400,
        wnae_pt_400_500,
        wnae_pt_500_inf,
        pt,
        wp,
    ):
        # High stat WP 17/11/2024 - 2018 but should work for all years
        if wp == 10:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 26.259)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 20.065)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 23.472)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 22.599)
               + (pt >= 500) * (wnae_pt_500_inf > 18.747)
            )
        elif wp == 20:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 25.156)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 18.284)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 20.383)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.941)
               + (pt >= 500) * (wnae_pt_500_inf > 16.370)
            )
        elif wp == 25:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.843)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.802)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 19.525)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.739)
               + (pt >= 500) * (wnae_pt_500_inf > 15.925)
            )
        elif wp == 30:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.594)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.440)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 18.776)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.558)
               + (pt >= 500) * (wnae_pt_500_inf > 15.468)
            )
        elif wp == 35:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.373)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.107)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 18.100)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.389)
               + (pt >= 500) * (wnae_pt_500_inf > 14.981)
            )
        elif wp == 40:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.166)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 16.781)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 17.493)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.218)
               + (pt >= 500) * (wnae_pt_500_inf > 14.526)
            )
        elif wp == 45:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 23.962)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 16.459)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 16.883)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.041)
               + (pt >= 500) * (wnae_pt_500_inf > 14.079)
            )
        else:
            print(f"Non valid working point: {wp}")
            exit(1)

        is_tagged = is_tagged.astype(int)
        #is_tagged = is_tagged.to_numpy().astype(int)
 
        return is_tagged

    j0_wnae_pt_0_200 = df["GoodJetsAK80_WNAEPt0To200Loss"]
    j1_wnae_pt_0_200 = df["GoodJetsAK81_WNAEPt0To200Loss"]
    j2_wnae_pt_0_200 = df["GoodJetsAK82_WNAEPt0To200Loss"]
    j3_wnae_pt_0_200 = df["GoodJetsAK83_WNAEPt0To200Loss"]

    j0_wnae_pt_200_300 = df["GoodJetsAK80_WNAEPt200To300Loss"]
    j1_wnae_pt_200_300 = df["GoodJetsAK81_WNAEPt200To300Loss"]
    j2_wnae_pt_200_300 = df["GoodJetsAK82_WNAEPt200To300Loss"]
    j3_wnae_pt_200_300 = df["GoodJetsAK83_WNAEPt200To300Loss"]

    j0_wnae_pt_300_400 = df["GoodJetsAK80_WNAEPt300To400Loss"]
    j1_wnae_pt_300_400 = df["GoodJetsAK81_WNAEPt300To400Loss"]
    j2_wnae_pt_300_400 = df["GoodJetsAK82_WNAEPt300To400Loss"]
    j3_wnae_pt_300_400 = df["GoodJetsAK83_WNAEPt300To400Loss"]

    j0_wnae_pt_400_500 = df["GoodJetsAK80_WNAEPt400To500Loss"]
    j1_wnae_pt_400_500 = df["GoodJetsAK81_WNAEPt400To500Loss"]
    j2_wnae_pt_400_500 = df["GoodJetsAK82_WNAEPt400To500Loss"]
    j3_wnae_pt_400_500 = df["GoodJetsAK83_WNAEPt400To500Loss"]

    j0_wnae_pt_500_inf = df["GoodJetsAK80_WNAEPt500ToInfLoss"]
    j1_wnae_pt_500_inf = df["GoodJetsAK81_WNAEPt500ToInfLoss"]
    j2_wnae_pt_500_inf = df["GoodJetsAK82_WNAEPt500ToInfLoss"]
    j3_wnae_pt_500_inf = df["GoodJetsAK83_WNAEPt500ToInfLoss"]

    j0_pt = df["GoodJetsAK80_pt"]
    j1_pt = df["GoodJetsAK81_pt"]
    j2_pt = df["GoodJetsAK82_pt"]
    j3_pt = df["GoodJetsAK83_pt"]

    svj0 = __is_tagged(
        j0_wnae_pt_0_200,
        j0_wnae_pt_200_300,
        j0_wnae_pt_300_400,
        j0_wnae_pt_400_500,
        j0_wnae_pt_500_inf,
        j0_pt,
        wp,
    )

    svj1 = __is_tagged(
        j1_wnae_pt_0_200,
        j1_wnae_pt_200_300,
        j1_wnae_pt_300_400,
        j1_wnae_pt_400_500,
        j1_wnae_pt_500_inf,
        j1_pt,
        wp,
    )

    svj2 = __is_tagged(
        j2_wnae_pt_0_200,
        j2_wnae_pt_200_300,
        j2_wnae_pt_300_400,
        j2_wnae_pt_400_500,
        j2_wnae_pt_500_inf,
        j2_pt,
        wp,
    )
    svj3 = __is_tagged(
        j3_wnae_pt_0_200,
        j3_wnae_pt_200_300,
        j3_wnae_pt_300_400,
        j3_wnae_pt_400_500,
        j3_wnae_pt_500_inf,
        j3_pt,
        wp,
    )
    
    n_svjs = svj0 + svj1 + svj2 + svj3

    return n_svjs



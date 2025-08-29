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

    ########################################
    ###   WNAE 2024_06_19 MET > 200 GeV  ###
    ########################################

    def __is_tagged(
        wnae_pt_0_200,
        wnae_pt_200_300,
        wnae_pt_300_400,
        wnae_pt_400_500,
        wnae_pt_500_inf,
        pt,
        wp,
    ):
        # "Low stat" WPs
        if wp == 10:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 83.15)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 29.60)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 14.34)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 39.98)
               + (pt >= 500) * (wnae_pt_500_inf > 10.91)
            )
        #elif wp == 20:
        #    is_tagged = (
        #       (pt < 200) * (wnae_pt_0_200 > 80.30)
        #       + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 27.97)
        #       + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 12.98)
        #       + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 35.78)
        #       + (pt >= 500) * (wnae_pt_500_inf > 9.79)
        #    )
        #elif wp == 30:
        #    is_tagged = (
        #       (pt < 200) * (wnae_pt_0_200 > 79.48)
        #       + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 26.92)
        #       + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 12.25)
        #       + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 32.77)
        #       + (pt >= 500) * (wnae_pt_500_inf > 9.45)
        #    )
        # "High stat" WPs
        elif wp == 20:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 79.554)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 28.145)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 13.008)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 35.670)
               + (pt >= 500) * (wnae_pt_500_inf > 9.779)
            )
        elif wp == 25:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 79.329)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 27.589)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 12.622)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 34.090)
               + (pt >= 500) * (wnae_pt_500_inf > 9.589)
            )
        elif wp == 30:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 79.177)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 26.956)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 12.293)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 32.714)
               + (pt >= 500) * (wnae_pt_500_inf > 9.433)
            )
        elif wp == 35:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 79.061)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 26.354)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 11.999)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 31.261)
               + (pt >= 500) * (wnae_pt_500_inf > 9.303)
            )
        elif wp == 40:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 78.972)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 25.722)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 11.732)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 29.917)
               + (pt >= 500) * (wnae_pt_500_inf > 9.182)
            )
        elif wp == 45:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 78.901)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 25.110)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 11.480)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 28.589)
               + (pt >= 500) * (wnae_pt_500_inf > 9.052)
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



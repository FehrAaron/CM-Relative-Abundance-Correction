# Correction for abundance changes (SW and "normal")


correct_abundance_changes = function(
    lip_data,
    trp_data,
    comparison,
    protein_id,
    grouping,
    diff,
    n_obs,
    std_error,
    p_adj_method = "BH",
    retain_columns = NULL,
    method = c("satterthwaite","no_df_approximation")){
  
  
  # Catch mistakes with selecting wrong method
  
  method = match.arg(method)
  
  
  # select required columns

    temp_lip_data = lip_data %>%
      dplyr::select(!!enquo(retain_columns), {{comparison}}, {{ protein_id }}, {{ grouping}}, {{ diff }}, {{ n_obs }}, {{ std_error }}) %>%
      dplyr::distinct()

    temp_trp_data = trp_data %>%
      dplyr::distinct({{comparison}}, {{ protein_id }}, {{ diff }}, {{ n_obs }}, {{ std_error }})

  # check for multiple entries per group:
    
    test = temp_lip_data %>%
      dplyr::distinct({{comparison}}, {{ protein_id }}, {{ grouping}})
    
    if (nrow(test) != nrow(temp_lip_data)){
      message("Your data frame contains dublicated values due to retained columns. This will affect the multiple testing correction.")
      }
    
    
  # define names for clearer code:

  se_pep =  rlang::sym(paste0(rlang::as_name(rlang::enquo(std_error)),"_pep"))

  se_prot = rlang::sym(paste0(rlang::as_name(rlang::enquo(std_error)),"_prot"))

  diff_pep = rlang::sym(paste0(rlang::as_name(rlang::enquo(diff)),"_pep"))

  diff_prot = rlang::sym(paste0(rlang::as_name(rlang::enquo(diff)),"_prot"))

  n_pep = rlang::sym(paste0(rlang::as_name(rlang::enquo(n_obs)),"_pep"))

  n_prot = rlang::sym(paste0(rlang::as_name(rlang::enquo(n_obs)),"_prot"))


# Join and correct data:

  # with satterthwaite:

  if (method == "satterthwaite"){

  corrected_data = dplyr::left_join(x = temp_lip_data,
                             y = temp_trp_data,
                             by = c(rlang::as_name(rlang::enquo(comparison)), rlang::as_name(rlang::enquo(protein_id))),
                             suffix = c("_pep", "_prot")) %>%
    dplyr::mutate(adj_diff = !!diff_pep - !!diff_prot) %>%
    dplyr::mutate(adj_std_error = sqrt( (!!se_pep)**2 + (!!se_prot)**2) ) %>%
    dplyr::mutate(numer = ((!!se_pep)**2 + (!!se_prot)**2)**2) %>%
    dplyr::mutate(denom = ( (!!se_pep)**4 / (!!n_pep-2) + (!!se_prot)**4/(!!n_prot-2))) %>%
    dplyr::mutate(df = numer/denom) %>%
    dplyr::mutate(tval = adj_diff / adj_std_error) %>%
    dplyr::mutate(pval = 2*stats::pt(abs(tval), df, lower.tail = FALSE))
  

  adjusted_data = dplyr::left_join(x = corrected_data,
                                   y = corrected_data %>%
                                     dplyr::filter(is.na(pval) == FALSE) %>%
                                     dplyr::group_by({{ comparison }}) %>%
                                     dplyr::mutate(adj_pval = p.adjust(pval, method = {{ p_adj_method }}))) %>%
    dplyr::select(- numer, - denom, - tval)

  return(adjusted_data)

  }

  # no df approximation
  
  
  if (method == "no_df_approximation"){

    corrected_data = dplyr::left_join(x = temp_lip_data,
                               y = temp_trp_data,
                               by = c(rlang::as_name(rlang::enquo(comparison)), rlang::as_name(rlang::enquo(protein_id))),
                               suffix = c("_pep", "_prot")) %>%

      dplyr::mutate(adj_diff = !!diff_pep - !!diff_prot) %>%
      dplyr::mutate(adj_std_error = sqrt( (!!se_pep)**2 + (!!se_prot)**2) ) %>%
      dplyr::mutate(df = !!n_pep - 2) %>%
      dplyr::mutate(tval = adj_diff / adj_std_error) %>%
      dplyr::mutate(pval = 2*stats::pt(abs(tval), df, lower.tail = FALSE))


    adjusted_data = dplyr::left_join(x = corrected_data,
                                     y = corrected_data %>%
                                       dplyr::filter(is.na(pval) == FALSE) %>%
                                       mutate(adj_pval = p.adjust(pval, method = {{ p_adj_method }}))) %>%
      dplyr::select(- tval)

    return(adjusted_data)

  }
  

}

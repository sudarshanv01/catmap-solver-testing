for d in coverages numbers_fix_xstar numbers_free_xstar coverages_ads_ads numbers_ads_ads_fix_xstar numbers_ads_ads_free_xstar; do
    echo $d
    cd $d
    bash cleanup.sh
    python mkm_job.py
    cd ..
    done

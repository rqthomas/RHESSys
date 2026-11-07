/******************************************************************************
 * rhessys_veg.c -- vegetation-backend selector (C side of the binding)       *
 *                                                                            *
 * Mirrors GLM3/src/glm_wqual.c (prime_wq). prime_veg() points the backend    *
 * function pointers at either the legacy (no-op here; big-leaf runs inline)  *
 * or the FATES Fortran BIND(C) entry points (rhessys_fates.F90).             *
 *                                                                            *
 * PHASE 0 (design doc 12): init/clean wired; set_bc/do are dummies until     *
 * Phases 1-3. Compile the FATES path with -DFATES and link libfates.         *
 ******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rhessys.h"
#include "rhessys_veg.h"

/*--------------------------------------------------------------------------*/
/* Legacy / dummy backend. The legacy big-leaf code is NOT a separate       */
/* module -- it runs inline in the daily cycle -- so these are no-ops and   */
/* exist only to give the pointers a safe default (off-mode, demog_flag==0).*/
/*--------------------------------------------------------------------------*/
static void legacy_veg_init(const char *pf, size_t *len, const char *inv, size_t *invlen,
                            const char *log, size_t *loglen, int *nlevsoil, int *nlevdecomp,
                            int *demog, int *nsites, int *is_restart)
{ (void)pf; (void)len; (void)inv; (void)invlen; (void)log; (void)loglen; (void)nlevsoil;
  (void)nlevdecomp; (void)demog; (void)nsites; (void)is_restart; }
static void dummy_veg_restart_write(const char *f, size_t *fl) { (void)f; (void)fl; }
static void dummy_veg_restart_read(const char *f, size_t *fl, struct veg_out *o)
{ (void)f; (void)fl; (void)o; }
static void dummy_veg_mgmt_init(const char *f, size_t *fl, int *n, int *p, int *z, int *h)
{ (void)f; (void)fl; (void)n; (void)p; (void)z; (void)h; }
static void dummy_veg_leaf_c(int *s, double *l) { (void)s; if (l) *l = 0.0; }
static void dummy_veg_inv_assign(int *n, int *i) { (void)n; (void)i; }
static void dummy_veg_set_ndemand(int *m, double *x) { (void)m; (void)x; }
static void dummy_veg_set_stomatal(int *m) { (void)m; }
static void dummy_veg_set_regen(int *m) { (void)m; }
static void dummy_veg_set_coage(int *m) { (void)m; }
static void dummy_veg_set_dayl(int *m) { (void)m; }
static void dummy_veg_set_disperse(int *m) { (void)m; }
static void dummy_veg_set_disperse_geom(int *n, double *x, double *y, double *a) { (void)n; (void)x; (void)y; (void)a; }
static void dummy_veg_set_l2fr_max(double *x) { (void)x; }
static void dummy_veg_set_min_n(double *x) { (void)x; }
static void dummy_veg_set_nvcmax(double *x) { (void)x; }
static void dummy_veg_set_nmine(double *x, int *m) { (void)x; (void)m; }
static void dummy_veg_set_nforage(double *x, int *m) { (void)x; (void)m; }
static void dummy_veg_set_nfun(double *k, int *c, int *m) { (void)k; (void)c; (void)m; }
static void dummy_veg_nfun_diag(int *s, double *v) { int k; (void)s; for (k = 0; k < 5; k++) v[k] = 0.0; }
static void dummy_veg_set_mort_disturb(int *s, double *f) { (void)s; (void)f; }
static void dummy_veg_set_browse(double *r, double *h, int *m) { (void)r; (void)h; (void)m; }
static void dummy_veg_set_suppress(double *g, double *d, int *n, double *r, int *l) { (void)g; (void)d; (void)n; (void)r; (void)l; }
static void dummy_veg_site_nneed(int *s, double *n) { (void)s; if (n) *n = 0.0; }
static void dummy_veg_litter_src(int *s, double *v) { int k; (void)s; for (k = 0; k < 10; k++) v[k] = 0.0; }
static void dummy_veg_set_leaf_export(int *s, double *f) { (void)s; (void)f; }
static void dummy_veg_leaf_export(int *s, double *v) { int k; (void)s; for (k = 0; k < 4; k++) v[k] = 0.0; }
static void dummy_veg_set_litter_direct(int *f) { (void)f; }
static void dummy_veg_litter_dcmp(int *s, double *v) { int k; (void)s; for (k = 0; k < 30; k++) v[k] = 0.0; }
static void dummy_veg_cwd_frac(double *v) { v[0] = 0.0; v[1] = 0.0; }
static void dummy_veg_init_site(int *site, int *patch_id, double *lat, struct veg_bc *bc)
{ (void)site; (void)patch_id; (void)lat; (void)bc; }
static void dummy_veg_coldstart(struct veg_out *out) { (void)out; }
static void dummy_veg_set_time(struct veg_time *t) { (void)t; }
static void dummy_veg_set_bc(int *site, struct veg_bc *bc, int *n_bad)
{ (void)site; (void)bc; if (n_bad) *n_bad = 0; }
static void dummy_veg_do(int *site, struct veg_bc *bc, struct veg_out *out)
{ (void)site; (void)bc; (void)out; }
static void dummy_veg_do_all(int *n, struct veg_bc *bcs, struct veg_out *outs)
{ (void)n; (void)bcs; (void)outs; }
static void dummy_veg_clean(void) {}
static void dummy_veg_dump(int *site, int *pid, int *ymd, const char *f, size_t *fl, int *first)
{ (void)site; (void)pid; (void)ymd; (void)f; (void)fl; (void)first; }

veg_init_t      p_veg_init      = legacy_veg_init;
veg_init_site_t p_veg_init_site = dummy_veg_init_site;
veg_coldstart_t p_veg_coldstart = dummy_veg_coldstart;
veg_set_time_t  p_veg_set_time  = dummy_veg_set_time;
veg_set_bc_t    p_veg_set_bc    = dummy_veg_set_bc;
veg_do_t        p_veg_do        = dummy_veg_do;
veg_do_all_t    p_veg_do_all    = dummy_veg_do_all;
veg_clean_t     p_veg_clean     = dummy_veg_clean;
veg_dump_t      p_veg_dump      = dummy_veg_dump;
veg_restart_write_t p_veg_restart_write = dummy_veg_restart_write;
veg_restart_read_t  p_veg_restart_read  = dummy_veg_restart_read;
veg_mgmt_init_t     p_veg_mgmt_init     = dummy_veg_mgmt_init;
veg_leaf_c_t        p_veg_leaf_c        = dummy_veg_leaf_c;
veg_inv_assign_t    p_veg_inv_assign    = dummy_veg_inv_assign;
veg_set_ndemand_t   p_veg_set_ndemand   = dummy_veg_set_ndemand;
veg_set_stomatal_t  p_veg_set_stomatal  = dummy_veg_set_stomatal;
veg_set_regen_t     p_veg_set_regen     = dummy_veg_set_regen;
veg_set_coage_t     p_veg_set_coage     = dummy_veg_set_coage;
veg_set_dayl_t      p_veg_set_dayl      = dummy_veg_set_dayl;
veg_set_disperse_t  p_veg_set_disperse  = dummy_veg_set_disperse;
veg_set_disperse_geom_t p_veg_set_disperse_geom = dummy_veg_set_disperse_geom;
veg_set_l2fr_max_t  p_veg_set_l2fr_max  = dummy_veg_set_l2fr_max;
veg_set_min_n_t     p_veg_set_min_n     = dummy_veg_set_min_n;
veg_set_nvcmax_t    p_veg_set_nvcmax    = dummy_veg_set_nvcmax;
veg_set_nmine_t     p_veg_set_nmine     = dummy_veg_set_nmine;
veg_set_nforage_t   p_veg_set_nforage   = dummy_veg_set_nforage;
veg_set_nfun_t      p_veg_set_nfun      = dummy_veg_set_nfun;
veg_nfun_diag_t     p_veg_nfun_diag     = dummy_veg_nfun_diag;
veg_set_mort_disturb_t p_veg_set_mort_disturb = dummy_veg_set_mort_disturb;
veg_set_browse_t   p_veg_set_browse   = dummy_veg_set_browse;
veg_set_suppress_t p_veg_set_suppress = dummy_veg_set_suppress;
veg_site_nneed_t    p_veg_site_nneed    = dummy_veg_site_nneed;
veg_litter_src_t    p_veg_litter_src    = dummy_veg_litter_src;
veg_set_leaf_export_t p_veg_set_leaf_export = dummy_veg_set_leaf_export;
veg_leaf_export_t   p_veg_leaf_export   = dummy_veg_leaf_export;
veg_set_litter_direct_t p_veg_set_litter_direct = dummy_veg_set_litter_direct;
veg_litter_dcmp_t   p_veg_litter_dcmp   = dummy_veg_litter_dcmp;
veg_cwd_frac_t      p_veg_cwd_frac      = dummy_veg_cwd_frac;

#ifdef FATES
/* FATES Fortran BIND(C) entry points (names match rhessys_veg.h _VEG_*).   */
void rhessys_veg_init(const char *paramfile, size_t *len, const char *invfile,
                      size_t *invlen, const char *logfile, size_t *loglen, int *nlevsoil,
                      int *nlevdecomp, int *demog_flag, int *nsites, int *is_restart);
void rhessys_veg_restart_write(const char *f, size_t *fl);
void rhessys_veg_restart_read(const char *f, size_t *fl, struct veg_out *o);
void rhessys_veg_mgmt_init(const char *f, size_t *fl, int *n, int *p, int *z, int *h);
void rhessys_veg_leaf_c(int *s, double *l);
void rhessys_veg_inv_assign(int *n, int *i);
void rhessys_veg_set_ndemand(int *m, double *x);
void rhessys_veg_set_stomatal(int *m);
void rhessys_veg_set_regen(int *m);
void rhessys_veg_set_coage(int *m);
void rhessys_veg_set_dayl(int *m);
void rhessys_veg_set_disperse(int *m);
void rhessys_veg_set_disperse_geom(int *n, double *x, double *y, double *a);
void rhessys_veg_set_l2fr_max(double *x);
void rhessys_veg_set_min_n(double *x);
void rhessys_veg_set_nvcmax(double *x);
void rhessys_veg_set_nmine(double *x, int *m);
void rhessys_veg_set_nforage(double *x, int *m);
void rhessys_veg_set_nfun(double *k, int *c, int *m);
void rhessys_veg_nfun_diag(int *s, double *v);
void rhessys_veg_set_mort_disturb(int *s, double *f);
void rhessys_veg_set_browse(double *r, double *h, int *m);
void rhessys_veg_set_suppress(double *g, double *d, int *n, double *r, int *l);
void rhessys_veg_site_nneed(int *s, double *n);
void rhessys_veg_litter_src(int *s, double *v);
void rhessys_veg_set_leaf_export(int *s, double *f);
void rhessys_veg_leaf_export(int *s, double *v);
void rhessys_veg_set_litter_direct(int *f);
void rhessys_veg_litter_dcmp(int *s, double *v);
void rhessys_veg_cwd_frac(double *v);
void rhessys_veg_coldstart(struct veg_out *out_all);
void rhessys_veg_init_site(int *site, int *patch_id, double *lat, struct veg_bc *bc);
void rhessys_veg_set_time(struct veg_time *t);
void rhessys_veg_set_bc(int *site, struct veg_bc *bc, int *n_bad);
void rhessys_veg_do(int *site, struct veg_bc *bc, struct veg_out *out);
void rhessys_veg_do_all(int *n, struct veg_bc *bcs, struct veg_out *outs);
void rhessys_veg_clean(void);
void rhessys_veg_dump_cohorts(int *site, int *pid, int *ymd, const char *f, size_t *fl, int *first);
#endif

/*--------------------------------------------------------------------------*/
int prime_veg(const char *which)
{
    if (which == NULL || strcmp(which, "legacy") == 0) {
        p_veg_init      = legacy_veg_init;
        p_veg_init_site = dummy_veg_init_site;
        p_veg_coldstart = dummy_veg_coldstart;
        p_veg_set_time  = dummy_veg_set_time;
        p_veg_set_bc    = dummy_veg_set_bc;
        p_veg_do        = dummy_veg_do;
        p_veg_do_all    = dummy_veg_do_all;
        p_veg_clean     = dummy_veg_clean;
        p_veg_dump      = dummy_veg_dump;
        p_veg_restart_write = dummy_veg_restart_write;
        p_veg_restart_read  = dummy_veg_restart_read;
        p_veg_mgmt_init     = dummy_veg_mgmt_init;
        p_veg_leaf_c        = dummy_veg_leaf_c;
        p_veg_inv_assign    = dummy_veg_inv_assign;
        p_veg_set_ndemand   = dummy_veg_set_ndemand;
        p_veg_set_stomatal  = dummy_veg_set_stomatal;
        p_veg_set_regen     = dummy_veg_set_regen;
        p_veg_set_coage     = dummy_veg_set_coage;
        p_veg_set_dayl      = dummy_veg_set_dayl;
        p_veg_set_disperse  = dummy_veg_set_disperse;
        p_veg_set_disperse_geom = dummy_veg_set_disperse_geom;
        p_veg_set_l2fr_max  = dummy_veg_set_l2fr_max;
        p_veg_set_min_n     = dummy_veg_set_min_n;
        p_veg_set_nvcmax    = dummy_veg_set_nvcmax;
        p_veg_set_nmine     = dummy_veg_set_nmine;
        p_veg_set_nforage   = dummy_veg_set_nforage;
        p_veg_set_nfun      = dummy_veg_set_nfun;
        p_veg_nfun_diag     = dummy_veg_nfun_diag;
        p_veg_set_mort_disturb = dummy_veg_set_mort_disturb;
        p_veg_set_browse   = dummy_veg_set_browse;
        p_veg_set_suppress = dummy_veg_set_suppress;
        p_veg_site_nneed    = dummy_veg_site_nneed;
        p_veg_litter_src    = dummy_veg_litter_src;
        p_veg_set_leaf_export = dummy_veg_set_leaf_export;
        p_veg_leaf_export   = dummy_veg_leaf_export;
        p_veg_set_litter_direct = dummy_veg_set_litter_direct;
        p_veg_litter_dcmp   = dummy_veg_litter_dcmp;
        p_veg_cwd_frac      = dummy_veg_cwd_frac;
        return 0;
    }

    if (strcmp(which, "fates") == 0) {
#ifdef FATES
        p_veg_init      = rhessys_veg_init;
        p_veg_init_site = rhessys_veg_init_site;
        p_veg_coldstart = rhessys_veg_coldstart;
        p_veg_set_time  = rhessys_veg_set_time;
        p_veg_set_bc    = rhessys_veg_set_bc;
        p_veg_do        = rhessys_veg_do;
        p_veg_do_all    = rhessys_veg_do_all;
        p_veg_clean     = rhessys_veg_clean;
        p_veg_dump      = rhessys_veg_dump_cohorts;
        p_veg_restart_write = rhessys_veg_restart_write;
        p_veg_restart_read  = rhessys_veg_restart_read;
        p_veg_mgmt_init     = rhessys_veg_mgmt_init;
        p_veg_leaf_c        = rhessys_veg_leaf_c;
        p_veg_inv_assign    = rhessys_veg_inv_assign;
        p_veg_set_ndemand   = rhessys_veg_set_ndemand;
        p_veg_set_stomatal  = rhessys_veg_set_stomatal;
        p_veg_set_regen     = rhessys_veg_set_regen;
        p_veg_set_coage     = rhessys_veg_set_coage;
        p_veg_set_dayl      = rhessys_veg_set_dayl;
        p_veg_set_disperse  = rhessys_veg_set_disperse;
        p_veg_set_disperse_geom = rhessys_veg_set_disperse_geom;
        p_veg_set_l2fr_max  = rhessys_veg_set_l2fr_max;
        p_veg_set_min_n     = rhessys_veg_set_min_n;
        p_veg_set_nvcmax    = rhessys_veg_set_nvcmax;
        p_veg_set_nmine     = rhessys_veg_set_nmine;
        p_veg_set_nforage   = rhessys_veg_set_nforage;
        p_veg_set_nfun      = rhessys_veg_set_nfun;
        p_veg_nfun_diag     = rhessys_veg_nfun_diag;
        p_veg_set_mort_disturb = rhessys_veg_set_mort_disturb;
        p_veg_set_browse   = rhessys_veg_set_browse;
        p_veg_set_suppress = rhessys_veg_set_suppress;
        p_veg_site_nneed    = rhessys_veg_site_nneed;
        p_veg_litter_src    = rhessys_veg_litter_src;
        p_veg_set_leaf_export = rhessys_veg_set_leaf_export;
        p_veg_leaf_export   = rhessys_veg_leaf_export;
        p_veg_set_litter_direct = rhessys_veg_set_litter_direct;
        p_veg_litter_dcmp   = rhessys_veg_litter_dcmp;
        p_veg_cwd_frac      = rhessys_veg_cwd_frac;
        fprintf(stderr, "RHESSys vegetation backend: FATES\n");
        return 0;
#else
        fprintf(stderr,
            "FATAL ERROR: -demog requested but RHESSys was built without FATES\n"
            "             (recompile with -DFATES and link libfates)\n");
        exit(EXIT_FAILURE);
#endif
    }

    fprintf(stderr, "FATAL ERROR: unknown vegetation backend '%s'\n", which);
    exit(EXIT_FAILURE);
}

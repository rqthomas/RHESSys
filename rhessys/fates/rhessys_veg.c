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
static void dummy_veg_set_ndemand(int *m) { (void)m; }
static void dummy_veg_site_nneed(int *s, double *n) { (void)s; if (n) *n = 0.0; }
static void dummy_veg_litter_src(int *s, double *v) { int k; (void)s; for (k = 0; k < 10; k++) v[k] = 0.0; }
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
veg_site_nneed_t    p_veg_site_nneed    = dummy_veg_site_nneed;
veg_litter_src_t    p_veg_litter_src    = dummy_veg_litter_src;
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
void rhessys_veg_set_ndemand(int *m);
void rhessys_veg_site_nneed(int *s, double *n);
void rhessys_veg_litter_src(int *s, double *v);
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
        p_veg_site_nneed    = dummy_veg_site_nneed;
        p_veg_litter_src    = dummy_veg_litter_src;
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
        p_veg_site_nneed    = rhessys_veg_site_nneed;
        p_veg_litter_src    = rhessys_veg_litter_src;
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

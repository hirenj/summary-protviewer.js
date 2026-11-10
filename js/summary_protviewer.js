import MASCP from 'mascp-jstools';
import 'sviewer/lite';
import {CondensedIupac, Sugar, SugarAwareLayoutFishEye, SVGRenderer} from 'glycan.js';

const ELEMENT_NAME = 'ccg-summary-protviewer';

const Iupac = CondensedIupac.IO;
const IupacSugar = Iupac(Sugar);

const GatorComponent = MASCP.GatorComponent;
const TrackRendererComponent = MASCP.TrackRendererComponent;
const Track = MASCP.Track;

const tmpl = document.createElement('template');

tmpl.innerHTML = `
<div id="summaries">
</div>
`
const summary_template = document.createElement('template');

summary_template.innerHTML = `
<input type="radio" name="summary"><label class="summary_group"><div class="sites"></div></label>
`

const entry_template = document.createElement('template');

entry_template.innerHTML = `
<div class="site_summary">
<ccg-sviewer-lite></ccg-sviewer-lite>
<div class="count"></div>
</div>
`


const style_tmpl = document.createElement('template');

style_tmpl.innerHTML = `
<style>
#summaries {
  --base-summary-size: 6em;
  --bottom-margin-size: 2px;
  --base-background: #eee;
  --selection-color: var(#f00,--selection-color);
  --button-default-background-color: var(#f00, --button-default-background-color);

  position: relative;
  width: 100%;
  /* Deliberately fixed, not content-driven: letting this grow/shrink with
     what's inside it (tried earlier) resizes the whole widget's layout,
     which can retrigger the very rendered/zoomdone/pandone events that
     rebuild this panel - a feedback loop that showed up as constant
     redrawing. An oversized composition is capped instead (see
     .site_summary's max-height) and just renders a bit smaller than its
     true per-residue scale, same graceful-degrade as a missing icon. */
  height: calc(var(--base-summary-size) + 1px + var(--bottom-margin-size));
  overflow: hidden;
}
.summary_group {
  display: flex;
  position: absolute;
  top: 0px;
  height: calc(100% - 1px - var(--bottom-margin-size));
  border-bottom: solid #999 var(--bottom-margin-size);

  flex-direction: column;
  justify-content: flex-end;
  align-items: center;
}

.summary_group.overflowed-content:after {
  display: flex;
  content: '+';
  justify-content: center;
  align-items: center;
  color: #fff;
  height: 1em;
  width: 1em;
  background: var(--button-default-background-color);
  box-shadow: 2px 2px 3px rgba(100,100,100,0.5);
  border-radius: 5px;
  position: absolute;
  top: 0px;
  left: 0px;
  cursor: pointer;
  z-index: 2;
}

.summary_group .sites {
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  align-items: flex-end;
  align-content: flex-end;
  justify-content: center;
  width: 100%;
  height: 100%;
}

.summary_group .sites {
  max-width: 100%;
}
.summary_group {
  z-index: 0;
}

input[type="radio"]:checked+.summary_group {
  z-index: 1;
  pointer-events: none;
}

input[type="radio"] {
  position: fixed;
  margin-left: -100vw;
}

input[type="radio"]:checked+.summary_group.overflowed-content .sites {
  width: 1000%;
  max-width: max-content;
}

input[type="radio"]:checked+.summary_group .sites .site_summary {
  /* Only the color changes here - box-sizing/border-width are set
     unconditionally on .site_summary itself so this never resizes it. */
  box-shadow: 2px 2px 3px rgba(100,100,100,0.5);
  border-color: var(--selection-color);
}

input[type="radio"]:checked+.summary_group .sites .site_summary .count {
  color: #222;
}

input[type="radio"]:checked+.summary_group.overflowed-content:after {
  display: none;
}

.site_summary {
  /* Defaults match the old fixed-size box exactly, so an entry looks
     the same as before until its composition's real size (see
     measure_natural_size/set_shape_vars) has been measured and applied.
     max-height caps it at the (fixed) panel budget - ccg-sviewer-lite
     scales its own content to fit whatever box it's given, so an
     oversized composition just renders a bit smaller than its true
     per-residue scale instead of being clipped or growing the panel. */
  --entry-height-fraction: 1;
  --entry-width-fraction: 1;
  height: calc( 0.5 * var(--base-summary-size) * var(--entry-height-fraction) );
  max-height: 100%;
  max-width: calc( 1.5 * var(--base-summary-size) );
  display: flex;
  flex-direction: row;
  background: var(--base-background);
  border-radius: 5px;
  /* Always present (transparent when not selected) and always border-box,
     so opening a group - which colors this border in, see the :checked
     rule below - never changes this box's actual pixel dimensions. It
     used to only gain a border (and switch box-sizing) once selected,
     shifting the content area by ~2px right at the moment of opening -
     enough for ccg-sviewer-lite to re-scale its internal SVG, which read
     as a flicker exactly when a group opened. */
  box-sizing: border-box;
  border: solid transparent 1px;
}
.site_summary .count {
  font-family: Helvetica, Arial, sans-serif;
  font-weight: bold;
  font-size: 9pt;
  color: #222;
  margin-right: 5px;
}
.site_summary ccg-sviewer-lite {
  --sugar-padding-side: 0.2;
  --sugar-padding-top: 0.2;
  height: 100%;
  width: calc( 0.5 * var(--base-summary-size) * var(--entry-width-fraction) );
}
</style>
`

const summarise_compositions = (sequence,groups) => {
  let max_peptide_composition = {};
  let site_compositions = {};
  for (let {index,ptm} of groups) {
    for (let {start,end,content,count} of ptm.filter( obj => obj instanceof Object)) {
      content = content.split('_').slice(1).join('_');
      max_peptide_composition[content] = Math.max(count,max_peptide_composition[content] || 0);
    }
    for (let content of ptm.filter( obj => typeof obj === "string")) {
      content = content.split('_').slice(1).join('_');
      if ( ! site_compositions[content]) {
        site_compositions[content] = {};
      }
      site_compositions[content][index] = true;
    }
  }
  return { peptides: max_peptide_composition, site_evidence: site_compositions };
}

const summarise_site_group = (group) => {
  let totals = group.map( site => site.ptm ).flat().map( comp => comp.split('_').slice(1).join('_'));
  let unique_entries = totals.filter( (o,i,a) => a.indexOf(o) == i );
  let tabled = Object.fromEntries(unique_entries.map( comp => {
    let matching = totals.filter( entry => entry == comp );
    return [comp, matching.length];
  }));

  return tabled;
};

const SEQUENCE_LOOKUP = {
  'gal(b1-3)galnac' : 'Gal(b1-3)GalNAc',
  'gal(b1-3)[glcnac(b1-6)]galnac' : 'Gal(b1-3)[GlcNAc(b1-6)]GalNAc',
  'glcnac(b1-3)[glcnac(b1-6)]galnac': 'GlcNAc(b1-3)[GlcNAc(b1-6)]GalNAc',
  'glcnac(b1-3)galnac' : 'GlcNAc(b1-3)GalNAc',
  'glcnac(b1-2)man(a1-3)[man(a1-6)]man(b1-4)glcnac(b1-4)glcnac' : 'GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
  'man(a1-2)man(a1-3)[man(a1-3)[man(a1-6)]man(a1-6)]man(b1-4)glcnac(b1-4)glcnac' : 'Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
  'man(a1-3)[man(a1-6)]man(b1-4)glcnac(b1-4)glcnac' : 'Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
  'galnac' : 'GalNAc',
  'glcnac' : 'GlcNAc',
  'man' : 'Man',
  'fuc' : 'Fuc',
  'xyl' : 'Xyl',
  'phospho' : 'P',
  // GAG-linker classes (see ccg-data-schemas/glycan-classes.json
  // GAG/CS/HS/KS, and the matching GAG_LINKER_STRUCTURES/KS_STRUCTURE in
  // ccg-shared-elements's static/msdata.packed.renderer.js) - without an
  // entry here, build_summary_elements() drops them from the summary
  // instead of just failing to render an icon for them.
  'galnac(b1-4)glca(b1-3)gal(b1-3)gal(b1-4)xyl' : 'GalNAc(b1-4)GlcA(b1-3)Gal(b1-3)Gal(b1-4)Xyl',
  'glcnac(b1-4)glca(b1-3)gal(b1-3)gal(b1-4)xyl' : 'GlcNAc(b1-4)GlcA(b1-3)Gal(b1-3)Gal(b1-4)Xyl',
  'hexnac(b1-4)glca(b1-3)gal(b1-3)gal(b1-4)xyl' : 'HexNAc(b1-4)GlcA(b1-3)Gal(b1-3)Gal(b1-4)Xyl',
  '{hso3(u1-6)glcnac(b1-3)[hso3(u1-6)]gal(b1-4)@y3a}n' : '{HSO3(u1-6)GlcNAc(b1-3)[HSO3(u1-6)]Gal(b1-4)@y3a}n'
}

// Each composition's natural (width,height) in glycan.js's own layout units
// - not just their ratio - measured once per unique composition (there are
// only ever as many as SEQUENCE_LOOKUP above has entries, so this is cheap
// in practice) by actually laying it out off-screen, the same way
// ensure_sugar_icon does in mascp-jstools's TrackRendererComponent.js - not
// guessed from the IUPAC string, since residue sizes and branch spacing
// aren't uniform enough for that to be reliable. Cache entries start as the
// in-flight Promise (so two entries needing the same composition share one
// render instead of racing) and get overwritten with the resolved
// {width,height} once it's ready, so later lookups are instant.
const natural_size_cache = new Map();

const measure_natural_size = (sequence) => {
  if (natural_size_cache.has(sequence)) {
    return Promise.resolve(natural_size_cache.get(sequence));
  }
  const promise = (async () => {
    // getBBox() (used by scaleToFit below) needs the element actually in the
    // rendered document to report real numbers in every browser - placed
    // off-screen rather than hidden via display/visibility, which breaks it
    // in some browsers (see the equivalent dance in ensure_sugar_icon).
    let holder = document.createElement('div');
    holder.style.position = 'absolute';
    holder.style.left = '-9999px';
    holder.style.top = '-9999px';
    document.body.appendChild(holder);
    try {
      let sug = new IupacSugar();
      sug.sequence = sequence;
      let renderer = new SVGRenderer(holder, SugarAwareLayoutFishEye);
      renderer.addSugar(sug);
      await renderer.refresh();
      renderer.scaleToFit({ side: 0.25, top: 0.5 });
      let viewBox = renderer.element.canvas.getAttribute('viewBox') || '';
      let [,,width,height] = viewBox.split(' ').map(parseFloat);
      return (width && height) ? {width,height} : {width: 1, height: 1};
    } catch (e) {
      return {width: 1, height: 1};
    } finally {
      holder.parentNode.removeChild(holder);
    }
  })();
  natural_size_cache.set(sequence, promise);
  promise.then( size => natural_size_cache.set(sequence, size) );
  return promise;
};

// Anchoring the scale on the smallest structure (a single residue) and
// letting bigger ones grow from there (tried first) can't give exact
// per-residue consistency while the panel itself stays a fixed height:
// anything tall enough eventually needs more height than the panel has and
// gets capped back down (by .site_summary's own max-height), which is
// exactly the same "smaller-than-true-scale box ->
// smaller residues" problem in a different place. The only way to
// guarantee every residue is genuinely the same size everywhere AND have
// everything fit is to anchor the other way: calibrate so the LARGEST
// known composition exactly fills the available space, and scale every
// other composition down from that same single factor. Smaller
// compositions render smaller than "today's default box" as a result -
// the explicit trade for exact consistency.
//
// Measured once, from every composition this component can ever show
// (SEQUENCE_LOOKUP is the complete, fixed set - nothing here is arbitrary
// user data), and cached for the component's whole lifetime.
let global_scale_promise = null;
const get_global_scale = () => {
  if ( ! global_scale_promise) {
    global_scale_promise = Promise.all(Object.values(SEQUENCE_LOOKUP).map(measure_natural_size))
      .then( sizes => 1 / Math.max(...sizes.map( ({height}) => height )) );
  }
  return global_scale_promise;
};

// A floor on width (and, symmetrically, height) specifically: a long
// linear chain with no branches (e.g. the GAG-linker classes - 5 residues
// in a straight line, see GAG_LINKER_STRUCTURES in msdata.packed.render
// er.js) is tall enough that true-to-scale width shrinks well past what
// ccg-sviewer-lite can actually contain cleanly, and the icon overflows its
// box instead of shrinking with it. Floor stops that at the cost of no
// longer being perfectly to-scale for the narrowest structures.
const MIN_FRACTION = 0.6;

const set_shape_vars = (el, size, scale) => {
  // No isConnected check here: this also runs synchronously from
  // apply_shape's cache-hit path, which fires from inside
  // build_summary_elements - *before* the cloned entry has been inserted
  // into the document at all. Checking isConnected here made every
  // cache-hit (i.e. every render after a composition's first) silently
  // skip setting the style, since the element was never connected yet at
  // the point this ran. Staleness only needs checking on the genuinely
  // async path (see apply_shape below), where the entry may have been
  // rebuilt or removed by the time the measurement resolves.
  let height_fraction = Math.max(MIN_FRACTION, size.height * scale);
  let width_fraction = Math.max(MIN_FRACTION, size.width * scale);
  el.style.setProperty('--entry-height-fraction', height_fraction.toFixed(3));
  el.style.setProperty('--entry-width-fraction', width_fraction.toFixed(3));
};

const apply_shape = (el, sequence) => {
  // Both go through cache/promise-sharing machinery of their own, so
  // asking for the scale every time costs nothing once it's warm.
  // Returned (not fire-and-forget) so callers - see build_summary_elements
  // below - can tell when every entry in a group has its real size, not
  // just the uniform default it starts out at.
  return Promise.all([measure_natural_size(sequence), get_global_scale()])
    .then( ([size, scale]) => {
      if ( ! el.isConnected) {
        // This entry's summary was rebuilt (or scrolled out and discarded)
        // before the off-screen measurement resolved - nothing to size.
        return;
      }
      set_shape_vars(el, size, scale);
    });
};

// Returns both the entry elements (append these, same as before) and
// `ready`, a Promise that resolves once every entry's real size has been
// applied - sizing is async (measure_natural_size), so .sites's layout
// right after appending these still reflects the uniform default size, not
// what it'll actually look like. Anything that needs to measure the real,
// final layout (see the overflow check in setSummaries) has to wait for
// this first.
const build_summary_elements = (summary) => {
  let pending = [];
  let entries = Object.entries(summary).map( ([content,count]) => {
    if ( ! SEQUENCE_LOOKUP[content]) {
      return;
    }
    let new_entry = entry_template.content.cloneNode(true);
    new_entry.querySelector('div.count').append(count);
    new_entry.querySelector('ccg-sviewer-lite').textContent = SEQUENCE_LOOKUP[content];
    pending.push(apply_shape(new_entry.querySelector('.site_summary'), SEQUENCE_LOOKUP[content]));
    return new_entry;
  }).filter( entry => entry);
  return { entries, ready: Promise.all(pending) };
}

const tally_compositions = function(sequence,peptides) {
  let positions = Object.values(peptides).flat().map( ({start,end}) => [start,end] ).flat();
  let all_contents = Object.values(peptides).flat().map( ({content}) => content.split('_').slice(1).join('_') ).filter( (o,i,a) => a.indexOf(o) == i );

  const min_pos = Math.min(...positions);
  const max_pos = Math.max(...positions);

  let all_densities = [];

  let densities_by_content = {};

  for (let {start,end,content,count} of Object.values(peptides).flat()) {
    let peptide_seq = sequence.substring(start-1, end);
    let possible_residues = (peptide_seq.match(/[STY]/g) || []).length;
    content = content.split('_').slice(1).join('_');
    let density = parseFloat((count/possible_residues).toFixed(2));
    all_densities.push({start,end,content,density});
  }
  for (let pos = min_pos ; pos <= max_pos; pos++) {
    let covering_peptides = all_densities.filter( ({start,end}) => pos >= start & pos <= end );
    for (let loop_content of all_contents) {
      if ( ! densities_by_content[loop_content]) {
        densities_by_content[loop_content] = new Array(sequence.length).fill(0);
      }
      densities_by_content[loop_content][pos] = Math.max(...[0].concat(covering_peptides.filter( ({content}) => content == loop_content ).map( ({density}) => density )));
    }
  }

  let results_pixels = {};
  let results_indexes = {};

  for (let loop_content of all_contents) {
    let densities = densities_by_content[loop_content];
    let splits = densities.map( (o,i,a) => { if (o > 0 && ! a[i-1]) return i - 1; if (!o && a[i-1]) return i - 1; }).filter( x => x );
    let splits_pixels = splits.reduce( (r,a,i) => {
      if (i % 2) {
        r[r.length - 1].push(this.screenPositionFor(a));
      } else {
        r.push([this.screenPositionFor(a)]);
      }
      return r;
    },[]);
    let splits_indexes = splits.reduce( (r,a,i) => {
      if (i % 2) {
        r[r.length - 1].push(a);
      } else {
        r.push([a]);
      }
      return r;
    },[]);

    results_pixels[loop_content] = splits_pixels;
    results_indexes[loop_content] = splits_indexes;
  }

  return {pixels: results_pixels, indexes: results_indexes};
};

const group_sites = function(sites) {
  let indexes = Object.keys(sites).map( aa => parseInt(aa));
  let pixels = indexes.map( aa => Math.floor(this.screenPositionFor(aa)) );
  let visible_positions = pixels.map( (o,i,a) => { 
    if (!Number.isNaN(o)) {
      return { pixel: o, index: indexes[i], ptm: sites[indexes[i]] }
    }
  }).filter( pos => pos );

  let pixels_diff = visible_positions.map( (o,i,a) => {
    if (i == 0) {
      return null;
    }
    return o.pixel - a[i-1].pixel;
  });

  const PIXEL_DISTANCE = 32;

  let splits = pixels_diff.map( (o,i) => { if (o > PIXEL_DISTANCE) return i; }).filter( x => x );
  let groups = splits.map( (split,idx,splits) => {
    let last_split = idx == 0 ? 0 : splits[idx-1];
    return visible_positions.slice(last_split,split);
  });
  groups.push(visible_positions.slice(splits[splits.length - 1]));
  return groups;
}

class SummaryViewer extends GatorComponent {
  connectedCallback() {
    super.connectedCallback();
    this.shadowRoot.querySelector('.widget_contents').prepend(tmpl.content.cloneNode(true));
    this.shadowRoot.append(style_tmpl.content.cloneNode(true));

    let composition_renderer = new TrackRendererComponent();
    const track_title = 'composition_summary';
    let new_track = new Track();
    new_track.setAttribute('name', track_title);
    new_track.setAttribute('fullname', 'Ambiguous');
    if ( ! this.querySelector(`*[name='${track_title}']`)) {
      this.appendChild(new_track);
    }
    composition_renderer.setAttribute('renderer', this.getAttribute('id'));
    composition_renderer.setAttribute('track', track_title);
    composition_renderer.script = function renderData(sequence,data){
      let results = [];
      let indexes_hexnac = new Array(sequence.length).fill(0);
      let indexes_hex = new Array(sequence.length).fill(0);

      for (let [composition,values] of Object.entries(data)) {
        let indexes = [];
        if (composition == 'hexnac' || composition == 'gal(b1-3)galnac') {
          indexes = indexes_hexnac;
        }
        if (composition == 'hex' || composition == 'man') {
          indexes = indexes_hex;
        }
        for (let [start,end] of values ) {
          while (start <= end) {
            indexes[start++ - 1] = 1;
          }
        }
      }
      if (indexes_hexnac.reduce((a, b) => a + b, 0) > 0) {
        results.push({ type: 'line',values: indexes_hexnac, options: { hide_axis: false, thickness: 2, color: '#ffd400', offset: -3 } });
      }
      if (indexes_hex.reduce((a, b) => a + b, 0) > 0) {
        results.push({ type: 'line',values: indexes_hex, options: { hide_axis: false, thickness: 2, color: '#00a651', offset: 3 } });
      }
      return results;
    };

    this.composition_renderer = composition_renderer;

    for (let renderer of this.ownerDocument.querySelectorAll(`ccg-trackrenderer[renderer="${this.getAttribute('id')}"]`)) {
      renderer.addEventListener('rendered', () => {
        this.updateSummaries();
      });
    }

    this.addEventListener('zoomdone', (ev) => {
      setTimeout(() => {
        this.updateSummaries();
      },50);
    });

    this.addEventListener('sequenceChange', (ev) => {
      this.updateSummaries();
    });

    this.addEventListener('pandone', (ev) => {
      this.updateSummaries();
    });
  }

  updateSummaries() {
    this.renderer.clearTrack(MASCP.getLayer('composition_summary'));
    this.setSummaries([]);
    for (let renderer of [...this.ownerDocument.querySelectorAll(`ccg-trackrenderer[renderer="${this.getAttribute('id')}"]`)]) {
      let visible = renderer.visible_items;
      if ( ! visible ) {
        continue;
      }
      let {ptms} = visible;
      if (! ptms ) {
        continue;
      }

      let peptides = Object.fromEntries(Object.entries(ptms).map( ([idx,ptm]) => {
        let filtered_ptms = ptm.filter( p => typeof p === 'object' );
        return [idx, filtered_ptms];
      }).filter( ([idx,entry]) => entry.length > 0 ));

      let sites = Object.fromEntries(Object.entries(ptms).map( ([idx,ptm]) => {
        let filtered_ptms = ptm.filter( p => typeof p === 'string' );
        return [idx, filtered_ptms];
      }).filter( ([idx,entry]) => entry.length > 0 ));

      let crowded_sites = group_sites.call(this,sites).filter( group => group.length > 0).filter( group => group.length > 1 || group[0].ptm.length > 1 );

      let { pixels: composition_summary_pixels, indexes: composition_summary_indexes } = tally_compositions.call(this,this.renderer.sequence,peptides);
      this.setSummaries(crowded_sites);
      this.composition_renderer.data = composition_summary_indexes;
    }
  }

  setSummaries(sitegroups) {
    let summaries = this.shadowRoot.querySelector('#summaries');
    while (summaries.firstChild) {
      summaries.removeChild(summaries.firstChild);
    }
    let radio_id=1;
    for (let group of sitegroups) {
      let site_summary = summarise_site_group(group);
      let new_entry = summary_template.content.cloneNode(true);
      let { entries: summary_elements, ready } = build_summary_elements(site_summary);
      new_entry.querySelector('label .sites').append(...summary_elements);
      new_entry.querySelector('label').setAttribute('for',`radio_${radio_id}`);
      new_entry.querySelector('input').setAttribute('id',`radio_${radio_id}`);
      radio_id++;
      summaries.appendChild(new_entry);
      new_entry = summaries.lastElementChild;
      let sites_el = new_entry.querySelector('.sites');
      if (sites_el.childNodes.length > 0) {
        // Measured directly (does the content need more height than is
        // actually available) once sizing has actually been applied, not
        // via an IntersectionObserver watching one particular child's
        // visibility (tried earlier, twice): that only reliably caught
        // overflow when the FIRST child specifically was the one being
        // clipped - with align-items:flex-end, a second, wrapped row can
        // overflow while the first (bottom) row stays fully visible, so a
        // group with genuinely more to show never got classified as
        // overflowing at all. It also reacted live to the group's own
        // open/close state, which changes .sites's width (see the
        // input:checked+.overflowed-content CSS rule) and so moved the
        // very child being observed - an open/close feedback loop once a
        // composition's size put the ratio right at the 0.1 threshold.
        //
        // Checking scrollHeight synchronously, right after appending,
        // doesn't work either: sizing is async (measure_natural_size), so
        // every entry is still at its uniform default size at that exact
        // moment and never looks like it overflows. `ready` (from
        // build_summary_elements) resolves once every entry in this group
        // has its real, final size applied - only then does .sites's
        // layout reflect what it'll actually look like.
        //
        // scrollHeight itself then turned out not to be reliable either -
        // WebKit/Safari in particular doesn't always report a flex-wrap
        // container's full wrapped content height through it. Checking
        // actual child geometry instead: align-items:flex-end anchors
        // content to the bottom, so overflow here specifically means a
        // wrapped child's box extends above .sites's own top edge
        // (clipped there by #summaries's overflow:hidden) - that's true
        // regardless of how any one browser computes scrollHeight.
        ready.then( () => {
          if ( ! new_entry.isConnected) {
            // This group was rebuilt (or removed) before sizing finished.
            return;
          }
          let top = sites_el.getBoundingClientRect().top;
          let overflowing = [...sites_el.children].some(
            child => child.getBoundingClientRect().top < top - 1
          );
          if (overflowing) {
            new_entry.classList.add('overflowed-content');
          }
        });
      }
      let added = summaries.lastElementChild;
      let min_pixel = Math.min(...group.map( gr => gr.pixel || 100000 ));
      let max_pixel = Math.max(...group.map( gr => gr.pixel || 0 ));
      if (min_pixel) {
        added.style.left = `${min_pixel}px`;
      }
      if (Math.abs(max_pixel - min_pixel) < 10) {
        added.style.left = `${Math.floor((max_pixel+min_pixel)/2)}px`;
        added.style.width = 'var(--base-summary-size)';//`${10}px`;        
        added.style.transform = 'translate(-50%)';
      } else if (max_pixel) {
        added.style.left = `${Math.floor((max_pixel+min_pixel)/2)}px`;
        added.style.width = `${max_pixel-min_pixel}px`;
        added.style.transform = 'translate(-50%)';
      }
    }
  }
}

customElements.define(ELEMENT_NAME,SummaryViewer);


export default SummaryViewer
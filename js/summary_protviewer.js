import MASCP from 'mascp-jstools';
import 'sviewer/js/sviewer/lite';

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
<x-sviewer-lite></x-sviewer-lite>
<div class="count"></div>
</div>
`


const style_tmpl = document.createElement('template');

style_tmpl.innerHTML = `
<style>
#summaries {
  --base-summary-size: 3em;
  --bottom-margin-size: 2px;
  --base-background: #eee;
  --selection-color: var(#f00,--selection-color);
  --button-default-background-color: var(#f00, --button-default-background-color);

  position: relative;
  width: 100%;
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
  box-shadow: 2px 2px 3px rgba(100,100,100,0.5);
  border: solid var(--selection-color) 1px;
  box-sizing: border-box;
}

input[type="radio"]:checked+.summary_group .sites .site_summary .count {
  color: #222;
}

input[type="radio"]:checked+.summary_group.overflowed-content:after {
  display: none;
}

.site_summary {
  height: calc( 0.5 * var(--base-summary-size) );
  max-width: var(--base-summary-size);
  display: flex;
  flex-direction: row;
  background: var(--base-background);
  border-radius: 5px;
}
.site_summary .count {
  font-family: Helvetica, Arial, sans-serif;
  font-weight: bold;
  font-size: 9pt;
  color: #222;
  margin-right: 5px;
}
.site_summary x-sviewer-lite {
  --sugar-padding-side: 0.2;
  --sugar-padding-top: 0.2;
  height: 100%;
  width: calc( 0.5 * var(--base-summary-size) );
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
  'man(a1-3)[man(a1-6)]man(b1-4)glcnac(b1-4)glcnac' : 'Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc',
  'galnac' : 'GalNAc',
  'glcnac' : 'GlcNAc',
  'man' : 'Man',
  'fuc' : 'Fuc',
  'xyl' : 'Xyl',
  'phospho' : 'P'
}

const build_summary_elements = (summary) => {
  return Object.entries(summary).map( ([content,count]) => {
    if ( ! SEQUENCE_LOOKUP[content]) {
      return;
    }
    let new_entry = entry_template.content.cloneNode(true);
    new_entry.querySelector('div.count').append(count);
    new_entry.querySelector('x-sviewer-lite').textContent = SEQUENCE_LOOKUP[content];
    return new_entry;
  }).filter( entry => entry);
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

    for (let renderer of this.ownerDocument.querySelectorAll(`x-trackrenderer[renderer="${this.getAttribute('id')}"]`)) {
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
    for (let renderer of [...this.ownerDocument.querySelectorAll(`x-trackrenderer[renderer="${this.getAttribute('id')}"]`)]) {
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
      let summary_elements = build_summary_elements(site_summary);
      new_entry.querySelector('label .sites').append(...summary_elements);
      new_entry.querySelector('label').setAttribute('for',`radio_${radio_id}`);
      new_entry.querySelector('input').setAttribute('id',`radio_${radio_id}`);
      radio_id++;
      summaries.appendChild(new_entry);
      new_entry = summaries.lastElementChild;
      if (new_entry.querySelector('.sites').childNodes.length > 0) {
        new IntersectionObserver(function([{intersectionRatio}]) {
          new_entry.classList.toggle('overflowed-content', intersectionRatio < 0.1);
        }, {root: new_entry}).observe(new_entry.querySelector('.sites').firstElementChild);
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

customElements.define('x-summary-protviewer',SummaryViewer);


export default SummaryViewer
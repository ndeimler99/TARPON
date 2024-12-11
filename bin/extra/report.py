"""Get a default report."""
from datetime import datetime
from typing import List, Optional, Type
import os
from dominate.tags import (
    a, button, code, div, h4, h2, html_tag, p, h1, section, attr, main, body, footer, head, header, style, title)

from ezcharts.components.params import ParamsTable
#from ezcharts.components.reports impor
from ezcharts.components.theme import (
    EPI2MELabsLogo, LAB_body_resources, LAB_head_resources)
from ezcharts.components.versions import VersionsTable
from ezcharts.layout.base import IClasses, Snippet
from ezcharts.layout.resource import Resource, ImageResource
from ezcharts.layout.snippets.banner import Banner
from ezcharts.layout.snippets.section import Section
from ezcharts.layout.util import cls, css, inline
from ezcharts.layout.base import IClasses, IStyles, Snippet

from typing import List, Type

from bokeh.embed import components
from dominate.util import raw

from ezcharts.components.ezchart import _BokehChart
from ezcharts.layout.base import Snippet
from ezcharts.layout.resource import (
    base_body_resources, base_head_resources, Resource)
from ezcharts.layout.snippets.document import DefaultBody, DefaultHead
from ezcharts.layout.util import write_report

class Report(Snippet):
    """A basic report."""

    TAG: str = 'html'

    def __init__(
        self,
        report_title,
        head_tag: Type[head] = DefaultHead,
        body_tag: Type[body] = DefaultBody,
        head_resources: List[Resource] = base_head_resources,
        body_resources: List[Resource] = base_body_resources
    ) -> None:
        """Create tag."""
        super().__init__(
            styles=None,
            classes=None)

        with self:
            # Adds the generic meta tags for us!
            self.head = head_tag()
            # Enables scroll spy globals for us!
            self.body = body_tag()

        with self.head:
            title(report_title)
            for resource in head_resources:
                resource()

        with self.body:
            self.header = header()
            self.main = main(style=css("background-color: #c8efdf; padding-top: 75px"))
            self.footer = footer()
            for resource in body_resources:
                resource()

    def write(self, path):
        """Write a report to file."""
        # check if the report contains `Bokeh` plots
        bokeh_charts = self.get_bokeh_charts()
        if bokeh_charts:
            # get the script + divs for all the plots; then place the div in the
            # corresponding `_BokehChart` div
            bokeh_script, bokeh_divs = components([x.plot._fig for x in bokeh_charts])
            for chart, bokeh_div in zip(bokeh_charts, bokeh_divs):
                with chart:
                    raw(bokeh_div)
            # add the script to the header
            with self.head:
                # make sure the plots fill out the enclosing div
                style(
                    """
                    .bk-Figure {
                        height: 100%;
                        width: 100%;
                    }
                    """
                )
                raw(bokeh_script)

        write_report(path, self)

    def get_bokeh_charts(self):
        """Return all children of the report that are of type `_BokehChart`."""
        bokeh_charts = []

        def _get_charts_in_children(s):
            if not hasattr(s, 'children') or not s.children:
                return
            for child in s.children:
                if isinstance(child, _BokehChart):
                    bokeh_charts.append(child)
                _get_charts_in_children(child)

        _get_charts_in_children(self)
        return bokeh_charts
    

class IBadgeClasses(IClasses):
    """Section html classes."""

    container: str = cls(
        "badge", "px-3", "py-2", "mb-2", "me-3", "badge")
    #container_bg: str = cls("bg-primary")
    container_bg: str = cls("py-2")

class IBadgeStyles(IStyles):
    """Section inline css styles."""

    container: str = css(
        "line-height: 20px", "border-radius: 6px",
        "font-size: 13px;")

    container2: str = css(
        "line-height: 20px", "border-radius: 6px", "margin-left: 15%",
        "font-size: 13px;")

    primary_badge: str = css("background-color: #214c46", "line-height: 20px", "border-radius: 6px",
        "font-size: 13px;")
    secondary_badge: str = css("background-color: #a8bdbf", "line-height: 20px", "border-radius: 6px",
        "font-size: 13px;")

class Badge(Snippet):
    """A styled span."""

    TAG = 'span'

    def __init__(
        self,
        title: str,
        bg_class=None,
        style=None,
        styles: IBadgeStyles = IBadgeStyles(),
        classes: IBadgeClasses = IBadgeClasses(),
    ) -> None:
        """Create styled badge."""
        bg = bg_class or classes.container_bg
        super().__init__(
            title,
            styles=styles,
            classes=classes,
            className=f"{classes.container} {bg}",
            style=style)

class ILabsAddendumClasses(IClasses):
    """Section html classes."""

    container: str = cls("py-5", "px-0", "border-top")
    inner: str = cls("container", "px-0")


class LabsAddendum(Snippet):
    """A styled footer component for use in a Report."""

    TAG = 'div'

    def __init__(
        self,
        workflow_name: str,
        workflow_version: str,
        classes: ILabsAddendumClasses = ILabsAddendumClasses(),
        use_defaults: bool = True
    ) -> None:
        """Create tag."""
        super().__init__(
            styles=None,
            classes=classes,
            className=classes.container)

        with self:
            self.container = div(className=classes.inner)

        if use_defaults:
            with self.container:
                h4('About this report', className="pb-3")
                p(
                    "This report was produced using the ",
                    code(f"{workflow_name}"),
                    f" nextflow workflow ({workflow_version}) with support of epi2me-labs/ezcharts packages."
                )
                p(
                    "Oxford Nanopore Technologies products are not "
                    "intended for use for health assessment or to "
                    "diagnose, treat, mitigate, cure or prevent any "
                    "disease or condition.")
                
                p(
                    "This workflow is not endorsed by, supported by, or associated with Oxford Nanopore Technologies."
                )


class ILabsNavigationClasses(IClasses):
    """Section html classes."""

    spacer: str = cls("d-flex")
    container: str = cls(
        "fixed-top", "d-flex", "align-items-center", "flex-wrap",
        "justify-content-center")
    inner: str = cls(
        "container", "px-0", "d-flex", "flex-wrap",
        "justify-content-center", "align-items-center", "py-2")
    logo: str = cls(
        "d-flex", "align-items-center", "pe-5", "mb-md-0",
        "me-md-auto", "text-decoration-none")
    dropdown_btn: str = cls("btn", "btn-primary", "dropdown-toggle")
    dropdown_menu: str = cls("dropdown-menu")
    dropdown_item_link: str = cls("dropdown-item")

class LabsNavigation(Snippet):
    """A styled nav component for use in a Report."""

    TAG = 'nav'

    def __init__(
        self,
        logo: Type[html_tag],
        groups: List[str],
        header_height: int = 75,
        classes: ILabsNavigationClasses = ILabsNavigationClasses()
    ) -> None:
        """Create tag."""
        spacer = div(
            className=classes.spacer,
            style=f"margin-top: {header_height}px;")

        super().__init__(
            styles=None,
            classes=classes,
            style=f"min-height: {header_height}px; background-color: #214c46",
            className=classes.container)

        spacer.add(self)
        with self:
            with div(className=self.classes.inner, styleName=IBannerStyles().container):
                with a(href="https://github.com/ndeimler99/TArPON", className=self.classes.logo):
                    h1("TArPON", className=IBannerStyles().inner, style=css("color: white"))

                button(
                    "Jump to section... ",
                    cls=self.classes.dropdown_btn,
                    type="button",
                    id="dropdownMenuButton",
                    data_bs_toggle="dropdown",
                    aria_haspopup="true",
                    aria_expanded="false",
                    style=css("background-color: #7bbcb6", "border-color:#7bbcb6"))

                ngroups = len(groups)
                with div(className=self.classes.dropdown_menu):
                    for count, group in enumerate(groups):
                        setattr(
                            self, group,
                            div(className='', __pretty=False))
                        if count != ngroups - 1:
                            div(cls="dropdown-divider")

    def add_link(
        self,
        group: str,
        link_title: str,
        link_href: str
    ) -> None:
        """Add a header nav link to the header links list."""
        group_list = getattr(self, group)
        with group_list:
            a(
                link_title,
                href=link_href,
                className=self.classes.dropdown_item_link)


class IBannerClasses(IClasses):
    """Section html classes."""

    #container: str = cls("px-0", "bg-dark", "labs-banner")

    container: str = cls("px-0")

    inner: str = cls(
        "container", "px-0", "py-2", "border-top", "text-white", "report-title")


class IBannerStyles(IStyles):
    """Section inline css styles."""

    container: str = css(
        "margin-bottom: -25px",
        "padding-bottom: 35px !important",
        #"background-color: #dc6a10")
        "background-color: #387f75")
    inner: str = css("border-color: rgba(255, 255, 255, 0.1) !important;")
        
class TarponLogo(div):
    """Labs logo element."""

    def __init__(self) -> None:
        """Create a div with an SVG logo inside."""
        print(os.getcwd())
        super().__init__(
            #inline(ImageResource('tarpon_logo.jpeg').data_file),
            inline("./tarpon_logo.jpeg"),
            tagname='div',
            style="width: 35px; height: 35px;",
            className="d-flex",
            alt="TArPON Logo")

class BasicReport(Report):
    """A basic labs-themed report."""

    def __init__(
        self,
        report_title,
        logo: Type[html_tag] = TarponLogo,
        head_resources: List[Resource] = LAB_head_resources,
        body_resources: List[Resource] = LAB_body_resources,
    ) -> None:
        """Create tag."""
        super().__init__(
            report_title=report_title,
           head_resources=head_resources,
           body_resources=body_resources
        )

        with self.header:
            self.nav = LabsNavigation(logo=logo, groups=['main', 'meta'])

        with self.main:
            self.intro_content = section(id="intro-content", role="region", style=css("background-color:#c8efdf"))
            self.main_content = section(id="main-content", role="region", style=css("background-color:#c8efdf"))


    def add_section(
        self,
        title: str,
        link: str,
        overflow: bool = False
    ) -> Section:
        """Add a section to the main_content region."""
        href = self.get_uid('Section')
        self.nav.add_link('main', link, f'#{href}')
        with self.main_content:
            return Section(href, title, overflow=overflow)


class LabsReport(BasicReport):
    """A basic labs-themed report for a workflow."""

    def __init__(
        self,
        report_title,
        workflow_name,
        workflow_params_path: str,
        workflow_versions_path: str,
        workflow_version: str,
        workflow_manifest_path: str,
        logo: Type[html_tag] = TarponLogo,
        head_resources: List[Resource] = LAB_head_resources,
        body_resources: List[Resource] = LAB_body_resources,
        created_date: Optional[str] = None
    ) -> None:
        """Create tag."""
        super().__init__(
            report_title=report_title,
            head_resources=head_resources,
            body_resources=body_resources)

        with self.header:
            self.nav.add_link('meta', 'Versions', '#versions')
            self.nav.add_link('meta', 'Parameters', '#parameters')
            self.nav.add_link('meta', 'Manifest', "#manifest")
            self.intro_content = section(id="intro-content", role="region")
            # with self.intro_content:
            #     self.banner = Banner(report_title, workflow_name)
            #     self.banner.add_badge("Research use only")
            #     if not created_date:
            #         created_date = datetime.today().strftime('%Y-%m-%d')
            #     self.banner.add_badge(created_date, bg_class="bg-secondary")
            #     self.banner.add_badge(workflow_version)
            with self.intro_content:
                with div(className=IBannerClasses().container, style=IBannerStyles().container):
                    h2(report_title, className=IBannerClasses().inner)
                    p(
                        f"Results generated through the {workflow_name} nextflow workflow",
                        className=IBannerClasses().inner
                    )
                    div.badges = div(className=IBadgeClasses().container, style=IBadgeStyles().container2)
                    with div.badges:
                        Badge("Research Use Only", style=IBadgeStyles().primary_badge)
                        Badge(datetime.today().strftime('%Y-%m-%d'), style=IBadgeStyles().secondary_badge)
                        Badge(workflow_version, style=IBadgeStyles().primary_badge)

        with self.main:
            self.meta_content = section(id="meta-content", role="region", style=css("background-color:#c8efdf"))
            with self.meta_content:
                with Section(
                    "versions",
                    'Software versions',
                    overflow=True
                ):
                    VersionsTable(workflow_versions_path)

                with Section(
                    "parameters",
                    'Workflow parameters',
                    overflow=True
                ):
                    ParamsTable(workflow_params_path)   

                with Section(
                    "manifest",
                    "Workflow Manifest",
                    overflow=True
                ):
                    ParamsTable(workflow_manifest_path)

        with self.footer:
            self.addendum = LabsAddendum(
                workflow_name=workflow_name, workflow_version=workflow_version
            )

    def add_badge(
        self,
        title: str,
        bg_class=None
    ) -> None:
        """Add a badge to the banner."""
        self.add_badge(title, bg_class=bg_class)
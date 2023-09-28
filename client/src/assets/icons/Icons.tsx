import BTree from './svg/BinaryTree.svg';
import Bone from './svg/bone.svg';
import Brain from './svg/brain.svg';
import Breast from './svg/breast.svg';
import Esophagus from './svg/esophagus.svg';
import Eye from './svg/eye.svg';
import Fly from './svg/fly.svg';
import Heart from './svg/heart.svg';
import HomoSapiens from './svg/homo_sapiens.svg';
import ImmuneSystem from './svg/immune_system.svg';
import Kidneys from './svg/kidneys.svg';
import Lung from './svg/lung.svg';
import MusMusculus from './svg/mus_musculus.svg';
import Ovary from './svg/ovary.svg';
import Rat from './svg/rat.svg';
import UserAnnotation from './svg/user_annotation.svg';
import CellTypeMarker from './svg/study_analysis.svg';
import ExpressionAnalysis from './svg/expression_analysis.svg';
import CoExpressionAnalysis from './svg/coexpression_analysis.svg';
import CompareAnnotations from './svg/annotation_comparison.svg';
import MSTeams from './svg/Microsoft_Teams.svg';

interface IIcon {
  className?: string;
  size?: string | number;
}

export function Icon({ className, size, icon, alt }: IIcon & { icon: string; alt: string }) {
  return <img src={icon} className={className} alt={alt} style={{ width: size || 'initial', height: size || 'initial' }} />;
}

export function BinaryTreeIcon({ className, size }: IIcon) {
  return <Icon icon={BTree} alt="binary tree" className={className} size={size} />;
}

export function BoneIcon({ className, size }: IIcon) {
  return <Icon icon={Bone} alt="bone" className={className} size={size} />;
}

export function BrainIcon({ className, size }: IIcon) {
  return <Icon icon={Brain} alt="brain" className={className} size={size} />;
}

export function BreastIcon({ className, size }: IIcon) {
  return <Icon icon={Breast} alt="breast" className={className} size={size} />;
}

export function EsophagusIcon({ className, size }: IIcon) {
  return <Icon icon={Esophagus} alt="esophagus" className={className} size={size} />;
}

export function EyeIcon({ className, size }: IIcon) {
  return <Icon icon={Eye} alt="eye" className={className} size={size} />;
}

export function FlyIcon({ className, size }: IIcon) {
  return <Icon icon={Fly} alt="fly" className={className} size={size} />;
}

export function HeartIcon({ className, size }: IIcon) {
  return <Icon icon={Heart} alt="heart" className={className} size={size} />;
}

export function HomoSapiensIcon({ className, size }: IIcon) {
  return <Icon icon={HomoSapiens} alt="homo sapiens" className={className} size={size} />;
}

export function ImmuneSystemIcon({ className, size }: IIcon) {
  return <Icon icon={ImmuneSystem} alt="immune system" className={className} size={size} />;
}

export function KidneysIcon({ className, size }: IIcon) {
  return <Icon icon={Kidneys} alt="kidneys" className={className} size={size} />;
}

export function LungIcon({ className, size }: IIcon) {
  return <Icon icon={Lung} alt="lung" className={className} size={size} />;
}

export function MouseIcon({ className, size }: IIcon) {
  return <Icon icon={MusMusculus} alt="mouse" className={className} size={size} />;
}

export function OvaryIcon({ className, size }: IIcon) {
  return <Icon icon={Ovary} alt="ovary" className={className} size={size} />;
}

export function RatIcon({ className, size }: IIcon) {
  return <Icon icon={Rat} alt="rat" className={className} size={size} />;
}

export function UserAnnotationIcon({ className, size }: IIcon) {
  return <Icon icon={UserAnnotation} alt="user annotation" className={className} size={size} />;
}

export function CellTypeMarkerIcon({ className, size }: IIcon) {
  return <Icon icon={CellTypeMarker} alt="cell type marker" className={className} size={size} />;
}

export function ExpressionAnalysisIcon({ className, size }: IIcon) {
  return <Icon icon={ExpressionAnalysis} alt="expression analysis" className={className} size={size} />;
}

export function CoExpressionAnalysisIcon({ className, size }: IIcon) {
  return <Icon icon={CoExpressionAnalysis} alt="coexpression analysis" className={className} size={size} />;
}

export function CompareAnnotationsIcon({ className, size }: IIcon) {
  return <Icon icon={CompareAnnotations} alt="compare annotations" className={className} size={size} />;
}

export function MSTeamsIcon({ className, size }: IIcon) {
  return <Icon icon={MSTeams} alt="microsoft teams" className={className} size={size} />;
}
